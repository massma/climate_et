#! ~/edward/bin/python

for i in range(5,7):
  print('working on fig %02d' % i)
  exec(open('fig%02d.py' % i).read())
exec(open('tables.py' % i).read())


