from src.processing.stream import Stream


class ForwardStream(Stream):
    def __init__(self, stream, hook=None, has_label=True):
        super().__init__(stream)
        self.stream = stream
        self.hook = hook
        self.has_label = has_label

    def __len__(self):
        return self.stream.__len__()

    def get(self, *args, **kwargs):
        element = self.stream.get(*args, **kwargs)

        if self.hook:
            if self.has_label:
                item = element[0]
            else:
                item = element

            self.hook(item)

        return element


class InverseMap(object):
    """Lookup map from feature to input (for traceability)."""

    def __init__(self):
        self.map = dict()
        self.ins = list()

    def input(self, stream, has_label=True):  # noqa: A003
        def hook(item):
            self.ins.append(item)

        return ForwardStream(stream, hook, has_label)

    def output(self, stream, has_label=True):
        def hook(item):
            assert len(self.ins) > 0
            key = item.tobytes()
            value = self.ins.pop(0)
            self.map[key] = value

        return ForwardStream(stream, hook, has_label)

    def find_input_for(self, output):
        return self.map.get(output.tobytes())


class NoOp(object):
    def input(self, stream, **kwargs):  # noqa: A003
        return stream

    def output(self, stream, **kwargs):
        return stream
