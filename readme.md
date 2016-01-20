# frobmask
tl;dr use lua 5.3!

just a little tool to find the frobenius number of a numerical semigroup. if your semigroup is <23,28,50> then

```
$ lua frobmask.lua 23 28 50
frobenius number is: 233
$ 
```

this program absolutely requires the integer features added in lua 5.3 so
make sure you are using a compatible interpreter
