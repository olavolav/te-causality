# TE-Causality - Example #1

This simple example should serve as an introduction on how to use TE-Causality.

The text file `signal.txt` contains to signals: The first one is simply white noise, drawn from a Gaussian distribution of mean zero and standard deviation of 1. The second on is also drawn (independently) from the same distribution, but then a fraction of the value of the first signal is added. In formulas:

    x_1(t) = eta_1(t)
    x_2(t) = eta_2(t) + 0.5*x_1(t-1)

As is done in the control file `control.txt`, all you need to do now is set a few basic parameters like `size` for the number of signals and `samples` for the number of samples available for each signal or node. `bins` is required by causality measures using binned estimates of probability distributions, namely `te-extended` and `mi`. And last but not least, we need to set `outputfile` and `outputparsfile` such that the output can be written to the right place in your file system.

Now we are good to go! You can compile any one of the causality measures included from the `transferentropy-sim/` directory. Type `rake -T` to get a list, or proceed to compile the Transfer Entropy measure like this:

    rake te-extended

Then you are ready to run the program, for instance by typing:

    ./te-extended examples/1/control.txt
