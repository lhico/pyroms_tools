import pkgutil
import importlib

for module_info in pkgutil.iter_modules(__path__):
    module = importlib.import_module(f"{__name__}.{module_info.name}")
    importlib.import_module(f"{__name__}.{module_info.name}")