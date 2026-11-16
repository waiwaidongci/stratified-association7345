from datetime import datetime
from sqlalchemy import String,ForeignKey,JSON,DateTime
from sqlalchemy.orm import Mapped,mapped_column
from src.config.database import Base
class SensitivityRun(Base):
    __tablename__='sensitivity_runs'
    id:Mapped[int]=mapped_column(primary_key=True)
    name:Mapped[str]=mapped_column(String(160),unique=True)
    base_run_id:Mapped[int]=mapped_column(ForeignKey('association_runs.id'))
    base_fingerprint:Mapped[str]=mapped_column(String(64),default='')
    swap_limit:Mapped[int]=mapped_column()
    status:Mapped[str]=mapped_column(String(20),default='queued')
    input_snapshot:Mapped[dict]=mapped_column(JSON,default=dict)
    result:Mapped[dict]=mapped_column(JSON,default=dict)
    error:Mapped[str]=mapped_column(String(700),default='')
    created_by:Mapped[int]=mapped_column(ForeignKey('users.id'))
    created_at:Mapped[datetime]=mapped_column(DateTime,default=datetime.utcnow)
