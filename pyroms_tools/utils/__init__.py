import importlib
import pkgutil

# Iterate over all modules in the current package
for module_info in pkgutil.iter_modules(__path__):
    module_name = module_info.name
    importlib.import_module(f'.{module_name}', package=__name__)