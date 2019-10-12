
Using the log files, I was able to construct a detailed investigation every 20,000 timesteps - so we have 100 pages to look through.

I've carried out a study of the reaction types and the lengths of the active and passive partners in each reaction. I think this is beggining to illumiate what is going on and tells us much about the cause of the floods. 

Open the pdf and I'll talk you through page 1 - the layout is the same for each page. The title of the plot has the timestep

There are six plots per page. The x-axis for each plot is the same, ranging from 0 to 100 opcodes

Remember, these are reactions *with no mutation* - but obviously the change from the seed parasite is due to mutation as well as the reaction dynamics shown here. The contribution of the raw copy mutation rate (i.e. when the `=` opcode is executed) to the observed dynamics will be a bit tricky to quantify...

The top two plots are barplots of the length distribution of the active and passive molecules for the timestep, you can see that by t=20000 we have a mode around length 64, which is the length of the original seed replicator. You can also see already that there are more short passive molecules. 

The remaining four plots take a little more explaining. I wanted to show which lengths are connected to which. First I grouped the reactions by type. Note that I've abandonned the distinction between 'self-self' and 'self-nonself' reactions for clarity. 

For each plot, there are two rows of dots. The top row shows the active molecules and the bottom row shows the passive molecules. Each reacting partner is connected with a line. I've coloured each reaction differently and made the colours transparent, to try and make it easier to see as many connections as possible. The size of each dot and width of each line is proportionate to the frequency in which a particular reaction is seen. Note this is on a log scale (plus a small constant) - so bigger dots are *much* more common than little dots. 

With this description in mind, we can begin to interpret what's going on. Lets start by looking at t=20000 on page 1. This is very early on in the evolution of the system. The first row shows the replicator reactions - those in which the passive partner is copied. You can see a big cluster of reactions around the modal length as we'd expect. 
