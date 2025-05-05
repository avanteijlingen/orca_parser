import sys, os
from .NWChemParse import *

if __name__ == "__main__":
    basename, ext = os.path.basename(sys.argv[1]).split(".")
    if ext == "nwout":
        oname = os.path.join(os.path.dirname(sys.argv[1]), f"{basename}.xyz")
        parser = NWChemParse(sys.argv[1])

        parser.parse()
        for frame in range(len(parser.coords)):
            mol = parser.makeASE(frame)
            mol.write(oname, append=(frame > 0))
        print(f"Written to: {oname}")
    else:
        print("Unknown file type")
