from pysdsl import _memory_monitor


class MemoryMonitor(object):

    def __init__(self, out_html=None, out_json=None):
        self.out_html = out_html
        self.out_json = out_json

    def __enter__(self):
        _memory_monitor.start()

    def __exit__(self, exc_type, exc_value, traceback):
        _memory_monitor.stop()

        if self.out_html is not None:
            _memory_monitor.report_html(self.out_html)

        if self.out_json is not None:
            _memory_monitor.report_json(self.out_json)
