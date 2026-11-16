import csv,io
from sqlalchemy import select
from fastapi import HTTPException
from src.models import AssociationRun,SensitivityRun
from src.config.queue import queue
from src.utils.audit import record
from src.utils.serializers import serialize
from .sensitivity_algorithm import analyze_sensitivity,plan_sensitivity
QUEUE='cenan:sensitivity'

def listing(db):
    rows=[]
    for x in db.scalars(select(SensitivityRun).order_by(SensitivityRun.id.desc())):
        r=serialize(x);r.pop('result');r.pop('input_snapshot');rows.append(r)
    return rows

def create(db,data,user):
    if db.scalar(select(SensitivityRun).where(SensitivityRun.name==data.name)):
        raise HTTPException(409,'敏感性复核名称重复')
    base=db.get(AssociationRun,data.run_id)
    if not base:raise HTTPException(404,'原研究不存在')
    if base.status!='completed' or not base.result:
        raise HTTPException(409,'仅已完成的原研究可做敏感性复核')
    snapshot=base.input_snapshot
    try:
        plan=plan_sensitivity({**snapshot,'alpha':base.alpha},data.swap_limit)
    except ValueError as exc:
        raise HTTPException(422,str(exc))
    frozen={'study':snapshot['study'],
            'tables':[{k:t[k] for k in ['code','name','a','b','c','d']} for t in plan['rows']],
            'alpha':base.alpha,'swap_limit':data.swap_limit,
            'source':{'base_run_id':base.id,'base_name':base.name,
                      'base_fingerprint':base.result.get('fingerprint',''),
                      'base_method_version':base.result.get('method_version','')}}
    row=SensitivityRun(name=data.name,base_run_id=base.id,
                       base_fingerprint=frozen['source']['base_fingerprint'],
                       swap_limit=data.swap_limit,input_snapshot=frozen,created_by=user.id)
    db.add(row);db.flush();record(db,user,'sensitivity.queued',row.id);db.commit()
    queue.rpush(QUEUE,str(row.id));return serialize(row)

def compute(db,id):
    row=db.get(SensitivityRun,id)
    if not row or row.status not in ['queued','running']:return
    base=db.get(AssociationRun,row.base_run_id)
    if not base or base.status!='completed' or not base.result:
        raise ValueError('原研究不再是已完成状态，拒绝复核')
    if base.result.get('fingerprint','')!=row.base_fingerprint:
        raise ValueError('原研究指纹与冻结时不一致，拒绝复核')
    row.status='running';db.commit()
    row.result=analyze_sensitivity(row.input_snapshot)
    _verify_baseline(row.result,base.result)
    if row.result['baseline']['alpha_tie']:
        # p exactly equals alpha: the frozen original run's stored flag is authoritative.
        inherited=bool(base.result['exact']['below_alpha'])
        row.result['baseline']['below_alpha']=inherited
        row.result['significance']['baseline_below_alpha']=inherited
    row.status='completed';row.error=''
    record(db,None,'sensitivity.completed',row.id);db.commit()
    return serialize(row)

def _verify_baseline(sensitivity,original):
    """Guard against method drift: the k=0 row must reproduce the original run."""
    s=sensitivity['baseline'];o=original
    if sensitivity['source'].get('base_fingerprint')!=o.get('fingerprint'):
        raise ValueError('敏感性快照记录的原研究指纹与冻结结果不一致')
    if s['observed_sum_a']!=o['observed_sum_a']:
        raise ValueError('敏感性基线Σa与原研究不一致')
    if abs(s['cmh']-o['cmh'])>1e-9*max(1.,abs(o['cmh'])):
        raise ValueError('敏感性基线CMH与原研究不一致')
    original_p=o['exact']['p_value']
    if original_p is not None and s['exact_p']['value'] is not None:
        if abs(s['exact_p']['value']-original_p)>1e-12*max(1.,abs(original_p)):
            raise ValueError('敏感性基线精确p与原研究不一致')
    if (s['adjusted_or']['status']==o['adjusted_or']['status']=='finite'
            and abs(s['adjusted_or']['value']-o['adjusted_or']['value'])>1e-12*max(1.,abs(o['adjusted_or']['value']))):
        raise ValueError('敏感性基线共同OR与原研究不一致')

def _ratio_text(x):
    return x['value'] if x['status']=='finite' else x['status']

def _cell(value):
    """Keep large integers exact in CSV instead of Excel float notation."""
    if isinstance(value,int) and not isinstance(value,bool) and abs(value)>=1_000_000_000_000:
        return f'{value:d}'
    return value

def csv_rows(run,note=None):
    """Render one frozen sensitivity result as CSV rows; never re-runs the enumeration."""
    r=run.result;b=r['baseline'];d=r['direction'];rows=[]
    def w(*values):rows.append([_cell(v) for v in values])
    w('sensitivity','method','fingerprint','base_run','base_fingerprint','alpha','swap_limit',
      'combination_count','combination_cap','support_size','boundary','tail_rule','numeric','baseline_note')
    w(run.name,r['method_version'],r['fingerprint'],run.base_run_id,run.base_fingerprint,r['alpha'],
      r['swap_limit'],r['capacity']['combination_count'],r['capacity']['combination_cap'],
      r['capacity']['support_size'],r['interpretation'],r['tail_rule'],r['numeric'],r['baseline_note'])
    if note is not None:
        w('note','title','interpretation','action','status','approved_by','review_comment')
        w('note',note.title,note.interpretation,note.action,note.status,note.approved_by,note.review_comment)
    w('baseline','direction','adjusted_or','cmh','asymptotic_p','exact_p','exact_p_log',
      'below_alpha','observed_sum_a')
    w('baseline',b['direction'],_ratio_text(b['adjusted_or']),b['cmh'],
      b['asymptotic_p']['display'],b['exact_p']['display'],b['exact_p']['log'],
      b['below_alpha'],b['observed_sum_a'])
    w('direction','higher','lower','equal','undefined','flip_count','changed_count',
      'examples_truncated')
    w('direction',d['counts']['higher'],d['counts']['lower'],d['counts']['equal'],
      d['counts']['undefined'],d['flip_count'],d['changed_count'],d['examples_truncated'])
    w('extreme','kind','tie_count','value','examples_truncated')
    for key in ['min','max']:
        w('extreme','exact_p_'+key,r['exact_p'][key]['tie_count'],
          r['exact_p'][key]['value']['display'],r['exact_p'][key]['examples_truncated'])
        w('extreme','cmh_'+key,r['cmh'][key]['tie_count'],repr(float(r['cmh'][key]['value'])),
          r['cmh'][key]['examples_truncated'])
    w('correction','kind','ordinal','observed_sum_a','direction','adjusted_or','cmh','exact_p',
      'exact_p_log','below_alpha','strata')
    groups=[('direction_flip',d['examples']),('exact_p_min',r['exact_p']['min']['examples']),
            ('exact_p_max',r['exact_p']['max']['examples']),('cmh_min',r['cmh']['min']['examples']),
            ('cmh_max',r['cmh']['max']['examples'])]
    for kind,items in groups:
        for x in items:
            w('correction',kind,x['ordinal'],x['observed_sum_a'],x['direction'],
              _ratio_text(x['adjusted_or']),x['cmh'],x['exact_p']['display'],
              x['exact_p']['log'],x['below_alpha'],
              ';'.join(f"{s['code']}:{s['k']:+d}" for s in x['strata']))
    w('per_layer','code','a','b','c','d','lower','upper','count')
    for p in r['capacity']['per_layer']:
        w('per_layer',p['code'],p['a'],p['b'],p['c'],p['d'],p['lower'],p['upper'],p['count'])
    return rows

def csv_body(run,note=None):
    out=io.StringIO();csv.writer(out).writerows(csv_rows(run,note));return out.getvalue()

def export(db,id,user):
    row=db.get(SensitivityRun,id)
    if not row:raise HTTPException(404,'敏感性复核不存在')
    if row.status!='completed' or not row.result:
        raise HTTPException(409,'敏感性复核尚未完成')
    body=csv_body(row)
    record(db,user,'sensitivity.exported',id);db.commit()
    return body
