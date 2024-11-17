#TODO: move to utils.utils
import yaml
import os

def load_yaml_with_includes(file_path, reference='default'):
    """
    Load a YAML file and any referenced YAML files recursively, merging them into a single dictionary.
    """
    with open(file_path, 'r') as file:
        config = yaml.safe_load(file)[reference]
    return resolve_references(config, os.path.dirname(file_path))


def resolve_references(config, base_path):
    """
    Recursively traverse the dictionary to load any referenced YAML files.
    """
    if isinstance(config, dict):
        for key, value in config.items():
            if isinstance(value, str) and value.endswith(('.yml', '.yaml')):
                # Resolve the full path to the referenced YAML file
                include_path = os.path.join(base_path, value) if not os.path.isabs(value) else value

                included_config = load_yaml_with_includes(include_path)
                config[key] = included_config
            else:
                # Recursively resolve nested dictionaries or lists
                config[key] = resolve_references(value, base_path)
    elif isinstance(config, list):
        return [resolve_references(item, base_path) for item in config]
    
    return config

#######################