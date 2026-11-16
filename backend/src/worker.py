import time,logging
from pathlib import Path
from src.config.database import SessionLocal
from src.config.queue import queue
from src.models import AssociationRun,SensitivityRun
from src.services.run_service import compute,QUEUE
from src.services.sensitivity_service import compute as compute_sensitivity, QUEUE as SENSITIVITY_QUEUE
logging.basicConfig(level=logging.INFO)
while True:
    Path('/tmp/heartbeat').touch()
    try:
        item=queue.blpop([QUEUE,SENSITIVITY_QUEUE],timeout=4)
        if not item:continue
        key,value=item
        model=AssociationRun if key==QUEUE else SensitivityRun
        with SessionLocal() as db:
            try:
                if key==QUEUE:compute(db,int(value))
                else:compute_sensitivity(db,int(value))
            except Exception as exc:
                db.rollback();run=db.get(model,int(value))
                if run:run.status='failed';run.error=str(exc)[:500];db.commit()
                logging.exception('association research failed')
    except Exception:logging.exception('worker loop');time.sleep(2)
