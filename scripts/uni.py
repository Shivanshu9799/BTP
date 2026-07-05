import timm
from timm.data import resolve_data_config
from timm.data.transforms_factory import create_transform

model = timm.create_model(
    "hf-hub:MahmoodLab/UNI",
    pretrained=True,
    init_values=1e-5,
    dynamic_img_size=True
)
transform = create_transform(**resolve_data_config(model.pretrained_cfg, model=model))
model.eval()
