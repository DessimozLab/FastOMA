


# https://docs.dask.org/en/stable/futures.html

futures = client.map(score, x_values)

best = -1
for future in as_completed(futures):
    y = future.result()
    if y > best:
        best = y

for future, result in as_completed(futures, with_results=True):

for batch in as_completed(futures, with_results=True).batches():
    for future, result in batch:
        ...

seq = as_completed(futures)

for future in seq:
    y = future.result()
    if condition(y):
        new_future = client.submit(...)
        seq.add(new_future)  # add back into the loop



# This is often used to signal stopping criteria or current parameters between clients.
var = Variable('stopping-criterion')
var.set(False)
var.get()
False


# If you want to share large pieces of information, then scatter the data first:

parameters = np.array(...)
future = client.scatter(parameters)
var.set(future)