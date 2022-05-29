def trim_zeros(array, loc='b', zero=0):
  assert loc in ['f', 'b', 'fb']
  start, end = 0, len(array)
  if 'f' in loc:
    while start < len(array) and array[start] == zero: start += 1
  if 'b' in loc:
    while end > 0 and array[end - 1] == zero: end -= 1
  return array[start:end]

def fill_zeros(array, n, loc='b', zero=0):
  assert loc in ['f', 'b']
  if loc == 'f':
    return [zero] * max(0, n - len(array)) + array
  if loc == 'b':
    return array + [zero] * max(0, n - len(array))