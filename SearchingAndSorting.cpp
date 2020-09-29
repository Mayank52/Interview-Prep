#include <iostream>
#include <vector>
#include <unordered_map>
#include <list>
#include <set>
#include <algorithm>

using namespace std;

// 169. Majority Element (Boyer-Moore Voting Algorithm)
/*
In case the majority element always exists, then sort the array and return arr[n/2]
Otherwise
Method 1- Brute Force
Method 2- Hashmap to keep freq of elements, and for each element check if their freq is >n/2: O(n),O(n)
Method 3- Use self balancing BST, similar to ordered map: O(nlogn), O(n), for normal BST O(n^2) time
Method 4- Sort the array, then keep track of count of elements: O(nlogn), O(1)
Method 5- Moore's Voting Algorithm: O(n), O(1)
*/
// Boyer-Moore Voting Algorithm
int majorityElement(vector<int> &nums)
{
    //find candidate
    int me = nums[0], count = 1;
    for (int i = 1; i < nums.size(); i++)
    {
        if (nums[i] == me)
            count++;
        else
            count--;

        if (count == 0)
        {
            me = nums[i];
            count = 1;
        }
    }

    //verify
    count = 0;
    for (int i = 0; i < nums.size(); i++)
        if (nums[i] == me)
            count++;

    if (count > nums.size() / 2)
        return me;

    return -1;
}
//Using Hashmap or BST
int majorityElement(vector<int> &nums)
{
    int n = nums.size();
    unordered_map<int, int> freq;
    for (int i = 0; i < n; i++)
    {
        freq[nums[i]]++;
        if (freq[nums[i]] > n / 2)
            return nums[i];
    }

    return -1;
}
//Incase it always exists
int majorityElement(vector<int> &nums)
{
    sort(nums.begin(), nums.end());
    return nums[nums.size() / 2];
}

// 229. Majority Element II
//Approach 1- Hashmap
vector<int> majorityElement2(vector<int> &nums)
{
    vector<int> res;
    int n = nums.size();
    unordered_map<int, int> freq;
    for (int i = 0; i < n; i++)
    {
        if (freq[nums[i]] != -1)
        {
            freq[nums[i]]++;
            if (freq[nums[i]] > n / 3)
            {
                res.push_back(nums[i]);
                freq[nums[i]] = -1;
            }
        }
    }

    return res;
}
//Approach 2- Modified Boyer Moore Voting Algorithm
vector<int> majorityElement2(vector<int> &nums)
{
    int n = nums.size();
    vector<int> res;
    int c1 = 0, c2 = 0, count1 = 0, count2 = 0;

    for (int i = 0; i < n; i++)
    {
        if (nums[i] == c1)
            count1++;
        else if (nums[i] == c2)
            count2++;
        else if (count1 == 0)
        {
            c1 = nums[i];
            count1 = 1;
        }
        else if (count2 == 0)
        {
            c2 = nums[i];
            count2 = 1;
        }
        else
        {
            count1--;
            count2--;
        }
    }

    count1 = 0;
    count2 = 0;
    for (int i = 0; i < n; i++)
    {
        if (nums[i] == c1)
            count1++;
        else if (nums[i] == c2)
            count2++;
    }

    if (count1 > n / 3)
        res.push_back(c1);
    if (count2 > n / 3)
        res.push_back(c2);

    return res;
}

// Searching in an array where adjacent differ by at most k
/*
A Simple Approach is to traverse the given array one by one and compare every element with given element ‘x’.
If matches, then return index.
The above solution can be Optimized using the fact that difference between all adjacent elements is at most k.
The idea is to start comparing from the leftmost element and find the difference between current array element and x.
Let this difference be ‘diff’. From the given property of array,
we always know that x must be at-least ‘diff/k’ away, so instead of searching one by one, we jump ‘diff/k’.
*/
int searchElement(vector<int> &nums, int k, int tar)
{
    int i = 0;
    while (i < nums.size())
    {
        // If x is found at index i
        if (nums[i] == tar)
            return i;

        // Jump the difference between current
        // array element and tar divided by k
        // We use max here to make sure that i
        // moves atleast one step ahead.
        i = i + max(1, abs((nums[i] - tar) / k));
    }

    //not found
    return -1;
}

// Find Missing And Repeating
/*
7 Methods-
https://www.geeksforgeeks.org/find-a-repeating-and-a-missing-number/

My method-
for all numbers mark arr[abs(arr[i])-1] -ve
if while marking -ve the target is already -ve, then that is repeating number
then for missing, iterate again, the index that has +ve element gives missing number
*/
int *findTwoElement(int *arr, int n)
{
    int missing, repeating;
    int *res = new int(2);
    for (int i = 0; i < n; i++)
    {
        if (arr[abs(arr[i]) - 1] < 0)
            repeating = abs(arr[i]);
        else
            arr[abs(arr[i]) - 1] *= -1;
    }

    for (int i = 0; i < n; i++)
    {
        if (arr[i] > 0)
        {
            missing = i + 1;
            break;
        }
    }

    res[0] = repeating;
    res[1] = missing;
    return res;
}

// Ceiling in a sorted array

// Find four elements that sum to a given value
/*
Approach 1-Recursion(Subsequence method) O(2^n)

Approach 2- O(NlogN), (N=n^2)
find all pair sums (n^2)
sort all pair sums, and use two pointer approach to find pairSum: O(NlogN)

Approach 3- O(n^2)
Make a hashmap of all pair sums
for elements in map, check if Sum-X exists in map, if it exists, then check if all 4 indexes are unique
*/
void fourElementSum(vector<int> &nums, int tar)
{
    int n = nums.size();
    unordered_map<int, vector<pair<int, int>>> mp; //{pair sum: list of indexes}
    set<vector<int>> res;

    //group all index pairs having same pair sum
    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
            mp[nums[i] + nums[j]].push_back({i, j});
    }

    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            //for each pair check give tar-pairSum exists in map
            int pairSum = nums[i] + nums[j];
            if (mp.find(tar - pairSum) != mp.end())
            {
                //if it exists, then pair it with all the indexes with that sum
                for (pair<int, int> ele : mp[tar - pairSum])
                {
                    int a = i;
                    int b = j;
                    int c = ele.first;
                    int d = ele.second;
                    //check that all indexes are unique
                    if (a != c && a != d && b != c && b != d)
                    {
                        //add them all in a vector and sort the vector
                        vector<int> smallAns;
                        smallAns.push_back(nums[a]);
                        smallAns.push_back(nums[b]);
                        smallAns.push_back(nums[c]);
                        smallAns.push_back(nums[d]);
                        sort(smallAns.begin(), smallAns.end());
                        //add it to a set so that duplicate ans are removed
                        res.insert(smallAns);
                    }
                }
            }
        }
    }

    if (res.size() == 0)
        cout << -1;

    for (auto ans : res)
    {
        for (auto ele : ans)
        {
            cout << ele << " ";
        }
        cout << "$";
    }
}


void solve()
{
}
int main()
{
    solve();
    return 0;
}