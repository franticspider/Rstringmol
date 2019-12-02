---
output:
  pdf_document: default
  html_document: default
---

# Introduction

I've carried out a study of the reaction types and the lengths of the active and passive partners in each reaction. I think this is beggining to illumiate what is going on and tells us much about the cause of the floods. 

Using the log files, I was able to construct a detailed investigation every 20,000 timesteps - so we have 100 pages to look through for each complete run. This sounds daunting, but if the pages are viewed as a pdf presentation, you can browse backward and forward through the runs easily and quickly.

# Reaction classification

This is given in detail in a separate document, but to summarise: 

- *replicator reactions* are reactions in which the passive partner is copied
- *parasitic reactions* are replicator reactions in which the replication is not reciprocated if the active/passive roles are reversed
- *reactions with no product* mean either the two reactants do not bind, or they bind and produce no new molecules
- *NonCatalytic reactions* are reactions in which one or both of the partner molecules are changed in the reaction 
- *Reactions with a different product*  leave the active and passive partners unchanged, but create new molecules that have a different sequence to the inputs. 

Remember, these are reactions *with no mutation* - but obviously the change from the seed parasite is due to mutation as well as the reaction dynamics shown here.

# Understanding the new plots

Open `smsp_out2.pdf` and I'll talk you through page 1 in detail -- the layout is the same for each page. Once I've done that, I'll point out how some of the events in this run can be interpreted using these plots.

The title of the plot has the run name and timestep. There are seven plots per page. From top to bottom we have:

- distribution of active species length
- distribution of passive species length
- Replicator reactions
- Parasite reactions
- Reactions with no product
- Non-catalytic reactions (macromutations)
- recations with a different product

The x-axis for each plot is the same, ranging from 0 to 100 opcodes

The top two plots are barplots of the length distribution of the active and passive molecules for the timestep, you can see that at t=20000 we have a mode around length 64, which is the length of the original seed replicator. You can also see already that there are more short passive molecules with length < 20. 

The remaining four plots take a little more explaining. I wanted to show which lengths are connected to which for each reaction class. First I grouped the reactions by type. Note that I've abandonned the distinction between 'self-self' and 'self-nonself' reactions here. 

For each plot, there are two rows of dots. The top row shows the active molecules and the bottom row shows the passive molecules. Each reacting partner is connected with a line.  The size of each dot and width of each line is proportionate to the frequency in which a particular reaction is seen. Note this is on a log scale (plus a small constant) - so bigger dots are *much* more common than little dots. I've coloured each reaction differently and made the colours transparent, to try and make it easier to see as many connections as possible.

With this description in mind, we can begin to interpret what's going on. Lets start by looking at t=20000 on page 1. This is very early on in the evolution of the system. The third plot row shows the replicator reactions - those in which the passive partner is copied. You can see a big cluster of reactions around the modal length as we'd expect. It's also clear that at this point, replicators are only binding to replicators of the same length - this is probably due to the spatial distribution of mutations in the early evolution of the system. 

The fourth plot row shows the parasites. The actual parasitic molecules are always the passive partners, so we can see the difference in length distributions between parasites and their hosts. The hosts are all the same length as the replicators in the plot above, but the parasite lengths are widely distributed. You can also see that some parasites are bound to some extremely long molecules beyond the range of the x axis

The fifth plot row show reactions with no product. It's perhaps surprising that there are so many of these. A lot are parasites bound to parasites - unable to copy each other. Some will be parasitic reactions in which the replicator has taken the passive role, so no product will be the result. Some will be species of functional replicators binding to other species of functional replicators. It'll take more work to tease this information out. 

The sixth row shows the noncatalytic reactions. There are fewer of these and again there is a cluster around the length of the functioning replicaors, indicating that these reactions occur when different species of replicator bind. 

The seventh row shows the reactions with a different product. These are rare, possibly due to the classification scheme I'm currently using. 

With this description in mind, let's now look at the dynamics of this run in detail. 

# Understanding the whole run

## General trends

If you look at the plot document as a presentation, you should be able to use the "page down" button to flip between timesteps quickly, which makes the plots approximate an animation of the entire run. Its worth doing this a few times to get familiar with the plots and how they change. The broad features are:

- Sudden jumps of the modal-length replicators to shorter/longer lengths
- Slow walk to shorter replicators on top of the sudden jumps

This summary shows that there is a mix of dynamics in the evolution of the system, but a more detailed examination of some of the features shows how different selection pressures are working at different times throughout the run. The youTube movies of the "birds eye view" of the arena for this run showed the following features: 

- Possible hypercycle at 12s,  t=290,000
- Coexistence of short and long replicator t = 300,000 to 420,000
- jumpers at t = 880,000
- flood at 48s, t = 1,136,000

## Hypercycles

Because we are focussing on length distributions, hypercycles are most evident in these plots where the lengths of the two partners differs. We aren't seeing this at t=280,000 or t=300,000 -- instead we see two distinct length distributions. It's possible that hypercycles exist for species of similar length, but they aren't clear in these plots. An example of candidate hypercycles can be seen at t=1,140,000, where there's clear mutual copying between molecules of lenth ~32 an those of length ~68, giving a characteristic 'X' shape in the plot of replicator reactions. 

## Short/long coexistence

Although hypercycles aren't clear, it is clear that two different size distributions coexist for a long time during this run. For many of the plots in this range, the longer species tolerate a large degree of parasitism, but the paraistes are commonly as long as their hosts. 

## Jumpers

These are shown as Non-Catalytic reactions, because one of the parent molecules is destroyed, resulting the the observed 'jumping' dynamics. The plot for t=880,000 captures this - a large number of molcules of length ~17 appear in this plot. 

## Flood

It's hard to pin a precise time on a flood event, so its worth comparing and scrolling between the state at t=1,260,000 and t= 1,460,000, either side of the notional flood event at t=1,360,000. Let's study the changes in each reaction type:

- *Replicators* There's a big increase in number, but there seems to be much more non-self replication, plus many hypercycles between different length replicators
- *Parasites*  are much reduced, particularly the short parasites (length <20) - this gives the momentum to the replicators, allowing the arena to fill up. Note also there's a small number of parasites on a hosts of length ~17 
- *No product* reactions also diminish, probably due to the high parasitic component of this class of reactions
- *NonCatalytic* reactions show a small increase, probably as a side effect of the increase in replicators. 

## After the flood

We see the following

- gradual shift of the modal replicator size from ~30 to ~20
- from t = 1,920,000, massive increase in long replicators at periodic size intervals - why?
- parastites stay at low number sun til arount t=1,700,000, then recover. 


# List of flood simulations

We now need to apply this analysis to the other runs:

*ranked by severity of flood*

- run 2 of 125x100, block seed at t = 1,100,000
- run 3 of 125x100, block seed at t = 1,640,000
- run 5 of 125x100, block seed at t =   580,000

*Other runs get dense, but don't really 'flood'*
