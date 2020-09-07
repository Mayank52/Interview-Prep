#include <iostream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

// 136. Single Number
int singleNumber(vector<int> &nums)
{
    int ans = nums[0];
    for (int i = 1; i < nums.size(); i++)
    {
        ans ^= nums[i];
    }
    return ans;
}

// 260. Single Number III
vector<int> singleNumber3(vector<int> &nums)
{
    int axorb = 0;
    for (int val : nums)
        axorb ^= val;

    int rsb = (axorb & -axorb);

    int a = 0, b = 0;
    for (int val : nums)
    {
        if ((rsb & val) == 0)
            a ^= val;
        else
            b ^= val;
    }

    return {a, b};
}

// 1342. Number of Steps to Reduce a Number to Zero
/*
For the binary representation from right to left(until we find the leftmost 1):
if we meet 0, result += 1 because we are doing divide;
if we meet 1, result += 2 because we first do "-1" then do a divide;
ony exception is the leftmost 1, we just do a "-1" and it becomse 0 already.

OR

To subtract 1 from an odd number means flipping a 1 to a 0.
To divide an even number by 2 means shifting right by one, dropping a 0.
After subtracting, you'll always get an even number
That means: (the number of steps required) = (the total number of bits in number) + (the number of set bits).
*/
int numberOfSteps(int num)
{
    if (num == 0)
        return 0;
    int count = 1;
    while (num > 1)
    {
        if ((num & 1) == 1)
            count += 2;
        else
            count += 1;
        num = num >> 1;
    }

    return count;
}

// 1404. Number of Steps to Reduce a Number in Binary Representation to One
/*
0+1=0
1+1=0, carry=1
*/
int numSteps(string s)
{
    int count = 0, c = 0;
    int i = s.size() - 1;
    while (i > 0)
    {

        if (s[i] == '0')
        {
            //if 0->even
            //c=1 means it becomes odd
            if (c == 1)
                count++;
        }
        else
        {
            //if 1->odd
            //c=1 means it becomes even, and c=0 means it remains odd and only in case of odd extra 1 is added to count
            if (c == 0)
                count++;
            c = 1;
        }
        count++;
        i--;
    }
    if (c == 1)
        count++;

    return count;
}

// 1442. Count Triplets That Can Form Two Arrays of Equal XOR
int countTriplets(vector<int> &arr)
{
    int triplets = 0;
    for (int i = 0; i < arr.size(); i++)
    {
        int xorSum = 0;
        for (int k = i; k < arr.size(); k++)
        {
            xorSum ^= arr[k];
            if (xorSum == 0)
                triplets += k - i;
        }
    }
    return triplets;
}

void solve()
{
}
int main()
{
    solve();
    return 0;
}