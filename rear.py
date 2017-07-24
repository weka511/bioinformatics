'''REAR 	Reversal Distance '''

def rear(s1,s2):
    print (s1,s2)
    def cycle(s,offset):
        return s[offset:]+s[:offset]
    def normalize(s,start=1):
        return cycle(s,s.index(start))
    def reverse(s,i0,i1):
        if i0<i1:
            return s[:i0] + s[i0:i1][::-1] + s[i1:]
        else:
            return reverse(cycle(s,i0),0,i1+len(s)-i0)
    def cyclically_equal(s1,s2):
        i=s2.index(s1[0])
        return s2==cycle(s1,-i)
    def generate_single_reversals(s):
        return [reverse(s,i0,i1) for i0 in range(len(s)) for i1 in range(len(s)) if i0!=i1]
    target=normalize(s1)

    history=set(str(normalize(s1)))
    candidates=[s2]
    for depth in range(25):
        print (depth)
        next_generation=[]
        for s in candidates:
            if target==normalize(s):
                return depth
            for rr in generate_single_reversals(s):
                key= str(normalize(rr))
                if not key in history:
                    next_generation.append(rr)
                    history.add(key)
        candidates=next_generation
    return None

def parse(line):
    return [int(c) for c in line.strip().split()]

#rear ([1,2,3,4,5],[1,2,3,4,5])
#rear ([1,2,3,4,5],[3,4,5,1,2])
#rear ([1,2,3,4,5],[3,4,5,2,1])
i=0

original=''

result=[]


with open (r'rear.txt') as f:
    for line in f:
        if i%3==0:
            original=parse(line)
        if i%3==1:
            result.append(rear(original,parse(line)))
        i+=1
    
print(result)
        