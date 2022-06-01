import functools
def multiplyLong(values):
    array = values.split(" ")
    for i, item in enumerate(array):
        array[i] = int(item)

    result = functools.reduce(lambda a,b : a*b, array)
    return result

value = input()
print(multiplyLong(value))