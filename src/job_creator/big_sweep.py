
parameter["shift"] = "shift1.txt"
#inputfile for generating the shift files
#~ gl = lambda x: range(0, shift.max_steps(x, x)+1)
#let the prog handle the spacing (uncomment line above and comment the two underneath)

parameter["spacing"] = 1.
gl = lambda x:[l for l in drange(0, x + parameter["spacing"] + 0.000001, parameter["spacing"])]
#set spacing manually

parameter["dirs"]  = [[l, gl(l)] for l in range(16, 17, 1)] 
#we only want one system with size 16
parameter["args"]   = [["L", "H"], "g"]
parameter["bash"]  = "-shift shift.txt -mult 1 " 
#-shift specifies the shift-file (since shift.txt is generated automatically leave it as is)
#-mult determines how many samples we want: mult == 1 <--> thermalization = 0.1M, samples = 1M
parameter["sq"]    = "-r 4h --mpp 500m"
#sharcnet specific commands, see/change lauch_programm
parameter["files"] = ["../../build/examples/sim"]
#where the executable is located
parameter["cmake"] = "-DUSE_S:STRING=2 -DUSE_GRID:STRING=3"
#USE_S is the renyi index and USE_GRID is the grid type (3=tri, 4=sqr, 6=hex)
#if you change this you need to recompile

