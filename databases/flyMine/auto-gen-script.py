#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("http://www.flymine.org/flymine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Interaction")

# The view specifies the output columns
query.add_view(
    "participant1.primaryIdentifier","participant2.primaryIdentifier",
    "details.experiment.publication.pubMedId",
    "details.experiment.interactionDetectionMethods.name"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Interaction.participant1.name", "ASC")

# You can edit the constraint values below
query.add_constraint("details.type", "=", "physical", code = "A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A")

## Anna added this at the bottom
out = open('flymine-query.txt','w')
out.write("#FBID1\tFBID2\tPubMedIDs\tExperiment\n")
for row in query.rows():
    out.write('%s\t%s\t%s\t%s\n' % (row["participant1.primaryIdentifier"], row["participant2.primaryIdentifier"], \
    row["details.experiment.publication.pubMedId"], row["details.experiment.interactionDetectionMethods.name"]))
out.close()
print 'Wrote to flymine-query.txt'