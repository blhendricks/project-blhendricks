## Astro 528 Class Project
## Background
Neutrinos are extremely difficult to detect due to their tiny mass and neutral charge. This also makes them fantastic messengers of distant sources since they're not deflected along the way. However, this makes them hard to detect directly. To get around this, in-ice detectors have been making use of radio waves that are produced via the Askaryan effect caused by neutrino interactions in the ice leading to a charge inbalance that travels faster than the local speed of ice in the medium. An important value for in-ice neutrino detectors is their effective volume, a quantitative measure of the sensitivity of the detector.

## Simulation
The effective volume of a detector can be calculated through simulation via Monte Carlo generation of neutrino effects along with detector information, ice properties, and the necessary neutrino physics. The first stage of this type of analysis is the generation of millions or more neutrino events to mimic experimental conditions. Various quantities have to be randomly sampled such as the minimum and maximum neutrino energy, arrival direction, and so on. With the high number of events necessary for accurate calculations, the Monte Carlo despite its strengths can be time-intensive. Thankfully, the method lends itself well to parallelization. The current simulation is not yet paralleized.

## Simulation Instructions
First clone the repository
```
$ git clone https://github.com/blhendricks/project-blhendricks.git
```

From inside the clone directory, type
```
$ chmod +x Event_Generator.jl
```
This gives you permission to run the file. Then, run
```
$ ./Event_Generator.jl
```
This will run the simulation with the default values. If you want to change the default values, use your preferred text editor to open the Event_Generator.jl file and change the final two lines to have the values you would like.
