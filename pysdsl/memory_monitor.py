from pysdsl import _memory_monitor_start, _memory_monitor_stop


class MemoryMonitor(object):

    def __init__(self):
        pass

    def __enter__(self):
        _memory_monitor_start()

    def __exit__(self, exc_type, exc_value, traceback):
        _memory_monitor_stop()
