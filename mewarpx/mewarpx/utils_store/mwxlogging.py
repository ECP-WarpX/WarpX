import logging

from mpi4py import MPI

comm = MPI.COMM_WORLD


class MEHandler(logging.StreamHandler):
    def __init__(self, stream):
        logging.StreamHandler.__init__(self, stream=stream)
        self.info_formatter = logging.Formatter(
            "%(message)s\n"
        )
        self.non_info_formatter = logging.Formatter(
            "%(levelname)s [%(name)s:%(lineno)s]: %(message)s\n"
        )

    def emit(self, record):
        stream = self.stream
        if record.levelname == 'INFO':
            stream.write(self.info_formatter.format(record))
        else:
            stream.write(self.non_info_formatter.format(record))


class MEFilter(logging.Filter):
    def filter(self, record):
        return comm.Get_rank() == 0
