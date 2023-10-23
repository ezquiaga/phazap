import os

_default_postprocessed_phase_dir = os.path.expanduser("~/.phazap")
if not os.path.exists(_default_postprocessed_phase_dir):
    try:
        os.makedirs(_default_postprocessed_phase_dir)
    except:
        # Fail silently
        _default_postprocessed_phase_dir = os.getcwd()