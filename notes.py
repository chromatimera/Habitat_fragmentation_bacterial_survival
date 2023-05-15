import operator

# Example 1: Use sorted() function to sort a list of lists
lists = [[93, 6], [72, 9], [35, 2]]
print("Sorted Lists based on index 0: % s" % (sorted(lists, key=operator.itemgetter(0))))
lists = [[2500, 'Spark'], [2200, 'Hadoop'], [3000, 'Python']]
print("Sorted Lists based on index 1: % s" % (sorted(lists, key=operator.itemgetter(1))))