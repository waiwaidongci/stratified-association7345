from pydantic import BaseModel,Field
class SensitivityInput(BaseModel):
    name:str=Field(min_length=3,max_length=160)
    run_id:int=Field(gt=0,strict=True)
    swap_limit:int=Field(default=1,ge=0,le=500,strict=True)
