--Copyright © 2015 Jon Allen <ylixir@gmail>
--This work is free. You can redistribute it and/or modify it under the
--terms of the Do What You Want To Public License,
--as published by Sam Hocevar. See the COPYING file for more details.

--this requires the bit twiddling features added in lua 5.3

--command line arguments should be the semigroup

--if the semigroup has a gcd, the program still plows on, looking
--for a frobenius number of the of the form n*gcd

--the basic data structure is an array of true false values.
--if the location in the array is in the semigroup then it
--is true. if it is not in the semigroup then it is false.

--this array is implemented by using a table of integers.
--every bit in every integer represents one true false value
--in our semigroup array.

--
--these are our variables
--

--the gcd of the semigroup
local base_divisor = tonumber(arg[1],10)
--the smallest and biggest number in the semigroup base
local base_min,base_max = base_divisor,base_divisor
--just arrays of integers for the bitmasks
local base_mask,group_mask = {},{}
--we also need a bitmask that will be the termination criteria
local term_mask
--the group_mask only keeps track of base_max numbers
--when we terminate, the frobenius number will be the
--number just before the start of our group_mask
local frobenius=-1
--the size of integers on this machine
local int_size=0

--
--these are our utility functions
--

--find out how many bits there are in an integer
local function getIntSize()
  local i=0
  repeat
    i=i+1
  until 0 == 1<<i
  return i
end

--if given an index of a number in the semigroup, then this
--will return the index of the integer in the table that
--contains the bit representing that number
local function int_index(i) return (i-1) // int_size + 1 end
--if given an index of a number in the semigroup, then this
--will return the index of the bit in the integer that contains it
local function bit_index(i) return (i-1) % int_size + 1 end

--seems like being able to find the gcd could be useful remember
--that ..and..or.. is the trinary operator in lua and so this is
--an absolute value embedded in the euclidean algorithm
local function gcd(x,y)
  return 0 == y and (x < 0 and -x or x) or gcd(y,x%y)
end

--this will generate a table that can be used to check for termination
--takes the smallest member of the semigroup basis and the semigroups
--greatest common divisor
local function make_terminator(min, div)
  local term = {}
  --fill in the 0s
  for i=1,int_index(min) do
    term[i] = term[i] or 0
  end
  --fill in the 1s
  for i=div,min,div do
    term[int_index(i)] = term[int_index(i)] | 1 << (i-1)%int_size
  end

  return term
end

--this will check and see if the routine is ready to terminate yet
local function check_terminator(group, term)
  for i,v in ipairs(term) do
    if v ~= v & group[i] then return false end
  end
  return true
end

--find the first non zero bit in a mask table
local function get_bit_offset(mask)
  for i,v in ipairs(mask) do
    if 0 ~= v then
      for j=1,int_size do
        if 0 ~= v & 1 << j-1 then return i,j end
      end
    end
  end
end
--note that this algorithm is awful, optimization
--should probably start here

--we need to know how many bits are in an integer
int_size=getIntSize()

--parse our command line, to get our semigroup
for i,v in ipairs(arg) do
  local v_int=tonumber(v,10)

  --a zero in our base is kind of meaningless...
  if 0 ~= v_int then
    --update our boundaries
    base_divisor = gcd(v_int,base_divisor)
    if v_int < base_min then base_min = v_int end
    if v_int > base_max then base_max = v_int end

    --build the base_mask, a bit array with a bit set for every
    --number in the semigroup base
    i = int_index(v_int)
    base_mask[i] = (base_mask[i] or 0) | 1 << bit_index(v_int)-1
  end
end

--I'm not interested in spending extra time branching in the actual
--algorithm. Let's fill in any nil spots in our table with zeros.
--we might as well get get our group generation started too
for i=1,int_index(base_max) do
  base_mask[i] = base_mask[i] or 0
  group_mask[i] = base_mask[i]
end

--we can edge case to keep from going past the end of the array,
--or we can just plug a zero in past the end so we don't worry
group_mask[#group_mask+1]=0

--the mask that we 'and' with our group mask to test for termination
term_mask = make_terminator(base_min, base_divisor)

--the meat, generate the semigroup in a window that goes from the
--conductor+1 to the conductor+base_mask until the bottom of the
--window contains our termination criteria
frobenius = -base_divisor
while false==check_terminator(group_mask, term_mask) do
  --get the first set bit
  i,j =  get_bit_offset(group_mask)

  --adjust our lower bound
  frobenius =  frobenius + (i-1)*int_size + j
  --shift our array down so the first  set bit drops off the bottom
  for k=i,#group_mask-1 do
    bottom = group_mask[k] >> j
    top = group_mask[k+1] << int_size - j
    group_mask[k-i+1] = bottom|top
  end
  --fill the end of our array with zeros
  for k=#group_mask-i+1,#group_mask-1 do group_mask[k]=0 end
  --add the new group elements in
  for k,v in ipairs(base_mask) do
    group_mask[k] = group_mask[k]|v
  end
end

print('frobenius number is: '..frobenius)

function check(num,list)
  if #list == 0 then
    if num == 0 then return {}
    else return false end
  else
    local new_list = {}
    for i=2,#list do new_list[i-1]=list[i] end
    for i=0,num//list[1] do
      local result = check(num-i*list[1],new_list)
      if result then
        result[list[1]]=i
        return result
      end
    end
    return false
  end
end

--[[ for to debug
local temp = {}
for i,v in ipairs(arg) do
  temp[i] = tonumber(v,10)
end
if check(frobenius,temp)==false then print (frobenius..' legit') end
for i=frobenius+base_divisor,frobenius+base_min,base_divisor do
  parts = check(i,temp)
  if false == parts then
    print(i..' not legit')
    break
  else
    local str = i..'='
    for j,v in pairs(parts) do
      str = str..j..'*'..v..'+'
    end
    print(str:sub(1,-2))
  end
end
]]
