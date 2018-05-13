source("Bool2Cont.r")

net = loadNetwork("Arabidopsis_thaliana_CC.net")
attr = getAttractors(net)
attr = Attr2df(attr)

Boolean_whole = ggplot(attr, aes(x = factor(Time), y = Variable, fill = factor(State)), colour = "black") + 
    geom_tile(width=0.9, height=0.9, size=3.0) +
    scale_fill_manual(values = c("#8DD35F", "#FF5555")) +
    #facet_grid(.~Attr, scales="free", space = "free") +  # Activate in case of multiple attractors
    theme_classic() +
    xlab("") +
    theme(axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
    labs(fill = "State") +
    theme(legend.title = element_text(face = "bold", size = 9.5)) +
    ggtitle("a) Boolean dynamic for whole A. thaliana cell cycle model") +
    theme(plot.title = element_text(size = 9, face = "bold"))

net = loadNetwork("backup.net")
attr = getAttractors(net)
attr = Attr2df(attr)

Boolean_simplified = ggplot(attr, aes(x = factor(Time), y = Variable, fill = factor(State)), colour = "black") + 
    geom_tile(width=0.9, height=0.9, size=3.0) +
    scale_fill_manual(values = c("#8DD35F", "#FF5555")) +
    # facet_grid(.~Attr, scales="free", space = "free") +  # Activate in case of multiple attractors
    theme_classic() +
    xlab("") +
    theme(axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()) +
    labs(fill = "State") +
    theme(legend.title = element_text(face = "bold", size = 9.5)) +
    ggtitle("c) Boolean dynamic for simplified A. thaliana cell cycle model") +
    theme(plot.title = element_text(size = 9, face = "bold"))


odeNet = BoolNet2ToContinuous("Arabidopsis_thaliana_CC.net")
times = seq(0, 35, 0.05)
o = ode(y = state, func = odeNet, parms = parameters, times = times)

o = as.data.frame(o)
o = melt(o, id = c("time"))
cont_whole = ggplot(subset(o, variable %in% c("KRP1", "CYCD31", "CYCB11")), aes(x = time, y = value, colour = variable)) + 
  geom_line() + 
  theme_linedraw() +
  ylab("Relative concentration") +
  xlab("Time (a.u.)") +
  theme(legend.title = element_text(face = "bold", size = 9.5)) +
  labs(colour = "Node", linetype = "Node") +
  ylim(0,1) + 
  scale_colour_brewer(palette="Set1") +
  ggtitle("b) Continuous dynamic for A. thaliana cell cycle model") +
  theme(plot.title = element_text(size = 9, face = "bold")) +
  theme(axis.title = element_text(size=7.5)) #+
  #facet_grid(variable~.)

odeNet = BoolNet2ToContinuous("backup.net")
times = seq(0, 35, 0.05)
o = ode(y = state, func = odeNet, parms = parameters, times = times)


o = as.data.frame(o)
o = melt(o, id = c("time"))
cont_simplified = ggplot(subset(o, variable %in% c("CYCD31", "CYCB11", "KRP1")), aes(x = time, y = value, colour = variable)) + 
  geom_line() + 
  theme_linedraw() +
  ylab("Relative concentration") +
  xlab("Time (a.u.)") +
  theme(legend.title = element_text(face = "bold", size = 9.5)) +
  labs(colour = "Node", linetype = "Node") +
  ylim(0,1) + scale_colour_brewer(palette="Set1") +
  ggtitle("d) Continuous dynamic for simplified A. thaliana cell cycle model") +
  theme(plot.title = element_text(size = 9, face = "bold")) +
  theme(axis.title = element_text(size=7.5)) 
  
  
grid.arrange(arrangeGrob(Boolean_whole + theme(legend.position="none"),
                         Boolean_simplified + theme(legend.position="bottom"),
                         nrow = 2, heights = c(1.75,1.25)), 
             arrangeGrob(cont_whole + theme(legend.position="none"),
                         cont_simplified + theme(legend.position="bottom"),
                         nrow = 2, heights = c(1, 1.25)),
                         ncol = 2)
