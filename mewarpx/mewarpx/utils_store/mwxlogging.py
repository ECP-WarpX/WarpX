import logging
from mpi4py import MPI

comm = MPI.COMM_WORLD


class MEHandler(logging.StreamHandler):
    def __init__(self, stream):
        logging.StreamHandler.__init__(self, stream=stream)
        self.formatter = logging.Formatter("%(name)s - %(levelname)s: %(message)s\n")

    def emit(self, record):
        stream = self.stream
        stream.write(self.formatter.format(record))


class MEFilter(logging.Filter):
    def filter(self, record):
        return comm.Get_rank() == 0
