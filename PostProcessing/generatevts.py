import paris
import os
import sys
extensions = ('.npz')

import pathlib

fname = sys.argv[1]
fout  = sys.argv[2]
if pathlib.Path(fout[:-4]+".vts").is_file() == False:
  paris.generate2DVTKNew(fname[:-4],fout[:-4])
if pathlib.Path("3D"+fout[:-4]+".vts").is_file() == False:
  paris.generate3DFrom2DVTKNew(fname[:-4],"3D"+fout[:-4])
if pathlib.Path("Ultra"+fout[:-4]+".vts").is_file() == False:
  paris.generate3DVTKUltra(fname[:-4],"Ultra"+fout[:-4])

