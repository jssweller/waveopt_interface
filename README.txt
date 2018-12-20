waveopt_interface.py

Genetic Algorithm optimizer for phase modulated wavefront. This program creates two win32 pipes for interfacing with SLM apparatus. For each generation each SLM mask is sent through output pipe to SLM apparatus and the output field metric for each mask is received from apparatus through input pipe. Masks and metrics are sent and received one at a time, sequentially. The metrics are then passed through a fitness function, ranking the masks which are sorted based on this rank.

The default fitness function simply takes the mean of the metric values received from the  pipe. The output field metric can be a single value or an array of values (e.g. pixel intensities in the region of interest).

Parameters can be adjusted with command line arguments. For a list of parameters that can be adjusted and their default values, run the following on the command line:

```bash
python waveopt_interface.py -h
```

For example, to run the optimizer for 500 generations, segment width of 16, and segment height of 12:

```bash
python waveopt_interface.py --gens 500 --segment_width 16 --segment_width 12
```
