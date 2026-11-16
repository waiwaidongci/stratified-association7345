from fastapi import APIRouter,Depends,HTTPException
from fastapi.responses import Response
from src.config.database import session
from src.constants.enums import Role
from src.middlewares.auth import require
from src.schemas.sensitivity import SensitivityInput
from src.controllers.sensitivity_controller import listing,create,export
from src.models import SensitivityRun,AssociationRun
from src.utils.serializers import serialize
router=APIRouter(prefix='/sensitivities')
@router.get('')
def index(db=Depends(session),u=Depends(require(Role.ADMIN,Role.ANALYST,Role.VIEWER))):return listing(db)
@router.get('/{id}')
def detail(id:int,db=Depends(session),u=Depends(require(Role.ADMIN,Role.ANALYST,Role.VIEWER))):
    row=db.get(SensitivityRun,id)
    if not row:raise HTTPException(404,'敏感性复核不存在')
    base=db.get(AssociationRun,row.base_run_id)
    return {**serialize(row),'base_name':base.name if base else ''}
@router.post('',status_code=202)
def add(data:SensitivityInput,db=Depends(session),u=Depends(require(Role.ADMIN,Role.ANALYST))):return create(db,data,u)
@router.get('/{id}/export')
def csv(id:int,db=Depends(session),u=Depends(require(Role.ADMIN,Role.ANALYST,Role.VIEWER))):
    return Response(export(db,id,u),media_type='text/csv; charset=utf-8',headers={'Content-Disposition':f'attachment; filename=cenan-sensitivity-{id}.csv'})
