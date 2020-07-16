import numpy 
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot 

# ========================== Loop Review =============================

# great for when you need to do one thing a lot of times 

# ====== While Loop 
# while <condition>:
    # do this

x = 0
while x < 10:
    print(x) 
    x = x + 1

# ====== For loop

# for each_entry in list_like:
    # do stuff. 

proton = ['up', 'up', 'down', 'gluon']
for item in proton:
    print(item)

# how many vowels are in the proton?
vowels = ['a', 'e', 'i', 'o', 'u']
vowels_counted = 0
for item in proton:
    for letter in item:
        if letter in vowels:
            vowels_counted += 1

# loop over list incl [0,1,...,9]. Print square of entry
for item in range(10):
    print(item**2)

# ========================== SOLUTIONS ================================

Events = numpy.loadtxt("Events.dat")

# Ex 1 - 2 minutes 

temp = Events[9][5]
print(temp)


# Ex 2 - 2 minutes 
temp = Events[1000][12]
print(temp)

# Ex 3 - 15 minutes 
# hist energy of highest pt Lep

# get list of energies
# so first we make two numpy arrays  which we'll use to store the energies 
number = len(Events)
energies_first  = numpy.zeros(number)
energies_second = numpy.zeros(number)

iterate = 0
while iterate < number:
    # then we loop over all the events and assign energies to our new arrays
    energies_first[iterate]  = Events[iterate][5]
    energies_second[iterate] = Events[iterate][11]

    iterate += 1

# make histograms 
bins = numpy.linspace(0,1000,101) #note! We use 101, since this is the number of edges and not bins! 

pyplot.hist(energies_first, bins)
pyplot.hist(energies_second, bins)
pyplot.xlabel("Energy [GeV]")
pyplot.ylabel("Count")
pyplot.show()

# Ex 4 - 15 minutes  

# Some counters we can use. Basically just integers. 
plus_plus   = 0
minus_minus = 0
plus_minus  = 0
minus_plus  = 0

iterate = 0
# we loop over all the events, and count them up! 
while iterate < number:
    event = Events[iterate]
    
    if event[6]>0:
        if event[12]>0:
            plus_plus += 1
        else:
            plus_minus += 1 
    else:
        if event[12]>0:
            minus_plus += 1
        else:
            minus_minus += 1

    iterate += 1

print("Both plus: {}".format(plus_plus))
print("Minus Plus: {}".format(minus_plus))
print("Plus Minus: {}".format(plus_minus))
print("Minus Minus: {}".format(minus_minus))

# Ex 5 - 15 minutes
# this one's a lot trickier! 
# we prepare some empty lists 
two_mu = []
two_el = []
mu_el  = []

iterate = 0
while iterate < number:
    # loop over the events
    # check what kind of event it is, then append the whole event to the relevant list 
    event = Events[iterate]
    if abs(event[7])==11 and abs(event[13])==11: # two electrons
        two_el.append(event)
    
    if abs(event[7])==13 and abs(event[13])==13: # two muons 
        two_mu.append(event)
    
    if (abs(event[7])==11 and abs(event[13])==13) or (abs(event[7])==13 and abs(event[13])==11):
        # electron and muon OR muon and electron
        mu_el.append(event)

    iterate += 1

print("Two Mu: {}".format(len(two_mu)))
print("Two El: {}".format(len(two_el)))
print("Combo: {}".format(len(mu_el)))

# Ex 6 - 3 minutes
# save the files
numpy.savetxt("two_mu.txt",two_mu)
numpy.savetxt("two_el.txt",two_el)
numpy.savetxt("mu_el.txt",mu_el)

### =========================== Alternate Solutions =====================================
print("Starting Alternates!")

# these are some other, weirder ways, we can solve the exercises

# Ex 3

# look up python list comprehension! It's fun 
energies  = numpy.transpose([[evt[5], evt[11]] for evt in Events])

pyplot.hist(energies[0], numpy.linspace(0,1000,100))
pyplot.hist(energies[1], numpy.linspace(0,1000,100))
pyplot.xlabel("Energy [GeV]")
pyplot.ylabel("Count")
pyplot.show()

# Ex 4

# Ex 5
two_mu = [event for event in Events if (abs(event[7])==13 and abs(event[13])==13) ]
two_el = [event for event in Events if (abs(event[7])==11 and abs(event[13])==11) ]
mu_el  = [event for event in Events if ((abs(event[7])==11 and abs(event[13])==13) or (abs(event[7])==13 and abs(event[13])==11))]

print("Two Mu: {}".format(len(two_mu)))
print("Two El: {}".format(len(two_el)))
print("Combo: {}".format(len(mu_el)))
