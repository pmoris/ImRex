import os
import time



class Model(object):
    def __init__(self, nameSuffix=""):
        self.name = "{}-{}".format(self.__class__.__name__, time.strftime("%Y%m%d-%H%M%S"))
        if nameSuffix:
            self.name += " - " + str(nameSuffix)

        self.name = self.name.strip()

    @property
    def baseName(self):
        return self.name

    def _buildModel(self):
        raise NotImplementedError

    def newInstance(self):
        return self._buildModel()

    def getName(self, iteration=None):
        name = self.name
        if iteration is not None:
            name += " ({})".format(str(iteration))
        return name

    def getOptimizer(self):
        raise NotImplementedError

    def getLoss(self):
        raise NotImplementedError

    @staticmethod
    def getCustomObjects():
        return {}
