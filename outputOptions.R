# These outputs all generate interactive "form" elements. They cannot be suspended
# because plots and tables depend on them

outputOptions(output, 'region', suspendWhenHidden=FALSE, priority=90)
outputOptions(output, 'cell.class', suspendWhenHidden=FALSE, priority=90)
outputOptions(output, 'cell.cluster', suspendWhenHidden=FALSE, priority=90)
outputOptions(output, 'cell.type', suspendWhenHidden=FALSE, priority=90)
outputOptions(output, 'current.cluster', suspendWhenHidden=FALSE, priority=100)
outputOptions(output, 'comparison.cluster', suspendWhenHidden=FALSE, priority=99)
outputOptions(output, 'current.subcluster', suspendWhenHidden=FALSE, priority=100)
outputOptions(output, 'comparison.subcluster', suspendWhenHidden=FALSE, priority=99)
outputOptions(output, 'dt.components.heading', suspendWhenHidden=FALSE, priority=99)

# These are the differential expression tables. They should render with a low
# priority so that the items above get computed first.
outputOptions(output, 'dt.cluster.markers', suspendWhenHidden=FALSE, priority=-100)
outputOptions(output, 'dt.subcluster.markers', suspendWhenHidden=FALSE, priority=-100)

