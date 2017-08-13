#################################################################################
# declare proxies in a separate file so I can interactively source the
# other .R files without errors.

dt.clusters.proxy = DT::dataTableProxy('dt.clusters')
dt.subclusters.proxy = DT::dataTableProxy('dt.subclusters')
