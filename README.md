# Comp-Physics-Project-A

This code should function automatically.
Simply run singlependulum.py or doublependulum.py

The simulation class also contains a module to generate animations like those
found in the ../animations/ folder. However, it was not able to run on the
windows computers in Blackett, so the relevant line: 

	sim.makeanimation(title=title)

is commented out. It may run on a linux machine, but in case not, the code runs without.

In addition, the singlependulum.py file has a variable 'find_thresholds'.
If this variable is set to True, the script will find the stepsize thresholds
for each method. Consequently, the code will take much longer to run. 
The default value is 'find_thresholds=false'.
