def whatprogram(mjd, programs, changes):
  index = 0
  
  if mjd >= changes[-1]:
    return  programs[-1]
  
  while mjd >= changes[index + 1]:
    index += 1

  return  programs[index]
