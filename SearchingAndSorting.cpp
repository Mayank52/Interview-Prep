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

// Find four elements that sum to a given value
/*
Approach 1-Recursion(Subsequence method) O(2^n)

Approach 2- (Best) O(NlogN), (N=n^2)
find all pair sums (n^2)
sort all pair sums, and use two pointer approach to find pairSum: O(NlogN)

Approach 3- O(n^2)
Make a hashmap of all pair sums
for elements in map, check if Sum-X exists in map, if it exists, then check if all 4 indexes are unique
to avoid duplicate quadruples, we use a set.
Now for an ordered_Set, complexity will be logn for insert, so overall complexity is not O(n^2)
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

// 18. 4Sum
vector<vector<int>> fourSum(vector<int> &nums, int tar)
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

    vector<vector<int>> ans;
    for (auto ele : res)
    {
        ans.push_back(ele);
    }

    return ans;
}

// Maximum sum such that no two elements are adjacent
/*
Use the recursion and apply dp
*/
long long maxSubsequenceSum(vector<int> &arr, int idx, vector<int> &dp)
{
    if (idx >= arr.size())
        return dp[idx] = 0;

    if (dp[idx] != -1)
        return dp[idx];

    //including the current element
    long long including = maxSubsequenceSum(arr, idx + 2, dp) + arr[idx];
    //excluding the current element
    long long excluding = maxSubsequenceSum(arr, idx + 1, dp);

    //return the max of inlcuding and excluding
    return dp[idx] = max(including, excluding);
}
long long maxSubsequenceSum_DP(vector<int> &arr)
{
    int n = arr.size();
    vector<int> dp(n + 2, -1);

    for (int idx = n + 1; idx >= 0; idx--)
    {
        if (idx >= arr.size())
        {
            dp[idx] = 0;
            continue;
        }

        //including the current element
        long long including = dp[idx + 2] + arr[idx];
        //excluding the current element
        long long excluding = dp[idx + 1];

        //return the max of inlcuding and excluding
        dp[idx] = max(including, excluding);
    }

    return dp[0];
}
long long maxSubsequenceSum(vector<int> &arr)
{
    int n = arr.size();
    vector<int> dp(n + 2, -1);
    return maxSubsequenceSum(arr, 0, dp);
}

// 349. Intersection of Two Arrays
/*
Approach 1- O(n+m) time, O(n) space
Use Hashmap to get count of elements in first array
then iterate over second array and check if it is present in map. 
If it is present, add it to result and remove it from map to avoid adding it again

Approach 2-O(1) space O(nlogn + mlogm) time
sort both arrays, and use two pointers to find the intersection
*/
vector<int> intersection(vector<int> &nums1, vector<int> &nums2)
{
    int n = nums1.size();
    int m = nums2.size();
    vector<int> res;
    sort(nums1.begin(), nums1.end());
    sort(nums2.begin(), nums2.end());

    int i = 0, j = 0;
    while (i < n && j < m)
    {
        //if it is equal , then add to result
        if (nums1[i] == nums2[j])
        {
            res.push_back(nums1[i]);
            //skip all other equal elements as it has been added already
            while (i < n - 1 && nums1[i] == nums1[i + 1])
                i++;
            while (j < m - 1 && nums2[j] == nums2[j + 1])
                j++;
            i++;
            j++;
        }
        else if (nums1[i] < nums2[j])
            i++;
        else
            j++;
    }

    return res;
}

// 350. Intersection of Two Arrays II
/*
Similar approaches as previous, just add the duplicates as well
*/
//Approach 1-Sorting (O(1) space)
vector<int> intersect(vector<int> &nums1, vector<int> &nums2)
{
    int n = nums1.size();
    int m = nums2.size();
    vector<int> res;

    sort(nums1.begin(), nums1.end());
    sort(nums2.begin(), nums2.end());

    int i = 0, j = 0;
    while (i < n && j < m)
    {
        if (nums1[i] == nums2[j])
        {
            res.push_back(nums1[i]);
            i++;
            j++;
        }
        else if (nums1[i] < nums2[j])
            i++;
        else
            j++;
    }

    return res;
}
//Approach 2- Hashmap (O(n) space)
vector<int> intersect(vector<int> &nums1, vector<int> &nums2)
{
    vector<int> res;
    unordered_map<int, int> freq;
    //count frequency of each element in nums1
    for (int ele : nums1)
        freq[ele]++;

    //for each element in nums2, if its count in map >0, add it to result, and decrease count by 1
    for (int ele : nums2)
    {
        if (freq[ele] > 0)
        {
            res.push_back(ele);
            freq[ele]--;
        }
    }

    return res;
}

// Find common elements in three sorted arrays (G4g test cases are wrong)
/*
Approach 1- Hashmap
Find intersection of first two, then find intersection of res and third array
Approach 2- Sorting
Using 3 pointers, similar to 2 pointers in 2 arrays
*/
vector<int> commonElements(int A[], int B[], int C[], int n1, int n2, int n3)
{
    vector<int> res;

    sort(A, A + n1);
    sort(B, B + n2);
    sort(C, C + n3);

    int i = 0, j = 0, k = 0;
    while (i < n1 && j < n2 && k < n3)
    {
        if (A[i] == B[j] && B[j] == C[k])
        {
            res.push_back(A[i]);
            while (i < n1 - 1 && A[i] == A[i + 1])
                i++;
            while (j < n2 - 1 && B[j] == B[j + 1])
                j++;
            while (k < n3 - 1 && C[k] == C[k + 1])
                k++;
            i++;
            j++;
            k++;
        }
        else if (A[i] < B[j])
            i++;
        else if (B[j] < C[k])
            j++;
        else
            k++;
    }

    return res;
}

// Count triplets with sum smaller than a given value
/*
Approach 1- Brute Force O(n^3)

Approach 2- O(nlogn + n^2)
Sort the array
Then use on for loop , for the first element
inside the for loop find the other two using 2 pointers approach
*/
long long countTriplets(long long arr[], int n, long long sum)
{
    sort(arr, arr + n);

    long long count = 0;

    for (int k = 0; k < n; k++)
    {
        int ele1 = arr[k];

        //use meet in the middle approach
        int i = k + 1, j = n - 1;
        while (i < j)
        {
            //sum is less so it is part of the ans
            //if for a fixed k, (i,j) gives sum(i,j,k)<target, then every element before j, will also give <sum with i,k
            //as in sorted array everything before j is smaller, so sum will keep decreasing for fixed i,k.
            if (arr[i] + arr[j] + arr[k] < sum)
            {
                //include all indexes from i+1 to in the ans, as they can all pair up with i,k
                count += j - i;
                i++;
            }
            //mySum>=target, so decrease j
            else
                j--;
        }
    }

    return count;
}

// Print all subarrays with 0 sum (count the subarrays)
int zeroSumSubarray(vector<int> &arr)
{
    unordered_map<long long, int> mp; //prefixSum:count
    long long currSum = 0;
    int count = 0;

    for (int i = 0; i < arr.size(); i++)
    {
        currSum += arr[i];
        if (currSum == 0)
            count++;
        if (mp.find(currSum) != mp.end())
            count += mp[currSum];
        mp[currSum]++;
    }

    return count;
}

// 462. Minimum Moves to Equal Array Elements II (Make all array elements equal with minimum cost)
/*
The quickest way to make them all equal is to make them all equal to the median
Median is the element in the middle of the sorted array
*/
int minMoves2(vector<int> &nums)
{
    int n = nums.size();

    //sort and get median
    sort(nums.begin(), nums.end());
    int median = nums[n / 2], count = 0;

    //find the count
    for (int i = 0; i < n; i++)
        count += abs(median - nums[i]);

    return count;
}


void solve()
{
}
int main()
{
    solve();
    return 0;
}