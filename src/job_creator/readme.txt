all this python script serve the purpose of seting up multiple (thousand) 
simulations in a conveniant way. 

Change the module-file (here big_sweep.cpp) to specify your sim (see comments in
that file) the module-file can be renamed.

When ready call "./generator" (perhaps you have to change the python version in 
the first line... "python-your-version ./generator" should always work) in the 
following way:

./generator module-file -cmd

for example:

./generator big_sweep.py -comp 

will compile the sim.cpp according to the specifications in the module file

./generator big_sweep.py -make 

will initialize the folder structure and copy the compiled sim in every folder
It also prepares the shift.txt files and bash_input.txt files.

./generator big_sweep.py -run

will call the command (that you need to change for your machine/cluster) 
on line 250 in function launch_program for all subfolders
that aren't done yet. It will not restart simulations that are done, but if you
call -run twice it will start the sim two times for each unfinished folder
(and you don't want that...). Make sure no sim runs for this module-file before
calling -run

./generator big_sweep.py -s (-a)

-s outputs the state of the sim. the -a shows the progress for all running sims.

(If you have multiple modules running, you can call ./generator -s (-a) to see
the progress on all modules)

./generator big_sweep.py -c -p

When the simulation is done, all subfolders will have a resutls.txt file.
-c collects them all into one colres.txt file and -p will plot according to 
the plotcommand in line 359 "bash_if("p", ...)"

./generator big_sweep.py -clean will delete all files belonging to the module
(so copy the plots/results somewhere else before) and deletes the sim executable
(you need to call -comp before -make works again)

if you have space issues you can call

./generator big_sweep.py -del

to delete just the sim executables in each subfolder. If you need them later again 
call -comp and -make again.
