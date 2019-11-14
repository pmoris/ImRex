from .stream import Stream


class ForwardStream(Stream):
    def __init__(self, stream, hook=None, hasLabel=True):
        super().__init__(stream)
        self.stream = stream
        self.hook = hook
        self.hasLabel = hasLabel

    def __len__(self):
        return self.stream.__len__()

    def get(self, *args, **kwargs):
        element = self.stream.get(*args, **kwargs)

        if self.hook:
            if self.hasLabel:
                item = element[0]
            else:
                item = element

            self.hook(item)

        return element


class InverseMap(object):
    def __init__(self):
        self.map = dict()
        self.ins = list()

    def input(self, stream, hasLabel=True):
        def hook(item):
            self.ins.append(item)

        return ForwardStream(stream, hook, hasLabel)

    def output(self, stream, hasLabel=True):
        def hook(item):
            assert len(self.ins) > 0
            key = item.tobytes()
            value = self.ins.pop(0)
            self.map[key] = value

        return ForwardStream(stream, hook, hasLabel)

    def findInputFor(self, output):
        return self.map.get(output.tobytes())


class NoOp(object):
    def input(self, stream, **kwargs):
        return stream

    def output(self, stream, **kwargs):
        return stream
