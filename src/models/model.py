import time


class Model(object):
    def __init__(self, name=None):
        if not name:
            self.name = "{}-{}".format(
                self.__class__.__name__, time.strftime("%Y-%m-%d-%H-%M-%S")
            )
        else:
            self.name = name
        self.name = self.name.strip()

    @property
    def base_name(self):
        return self.name

    def _build_model(self):
        raise NotImplementedError

    def new_instance(self):
        return self._build_model()

    def get_name(self, iteration=None):
        name = self.name
        if iteration is not None:
            name += " ({})".format(str(iteration))
        return name

    def get_optimizer(self):
        raise NotImplementedError

    def get_loss(self):
        raise NotImplementedError

    @staticmethod
    def get_custom_objects():
        return {}
