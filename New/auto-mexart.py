from sys import argv
from AutoCIP import PyCIP

if __name__  == "__main__":
    if len(argv) > 1 and all(tuple(map(lambda x: x.endswith(".dat")))):
        Files = argv[1:]
        CIP_Analysis = PyCIP(Files)