#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <queue>

using namespace std;

// 287. Find the Duplicate Number
int findDuplicate(vector<int> &nums)
{
    int slow = nums[0];
    int fast = nums[0];

    do
    {
        slow = nums[slow];
        fast = nums[nums[fast]];
    } while (slow != fast);

    slow = nums[0];
    while (slow != fast)
    {
        slow = nums[slow];
        fast = nums[fast];
    }

    return slow;
}

//KCON (Codechef)
long long maxSubarraySum(int arr[], int n)
{
    long long maxSoFar = arr[0];
    long long currMax = 0;

    for (int i = 0; i < n; i++)
    {
        currMax += arr[i];
        maxSoFar = max(maxSoFar, currMax);
        if (currMax < 0)
            currMax = 0;
    }

    return maxSoFar;
}
void KCON()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int n, k;
        cin >> n >> k;
        int a[n], b[n * 2];
        for (long long i = 0; i < n; i++)
        {
            cin >> a[i];
            b[i] = b[i + n] = a[i];
        }

        //if k is 1 then simply print max subarray sum
        if (k == 1)
        {
            cout << maxSubarraySum(a, n) << endl;
            continue;
        }

        //calculate maxSubarray sum for A*2 i.e A put 2 times
        long long maxSum = maxSubarraySum(b, n * 2);
        long long maxPref = -1e9, maxSuff = -1e9, currPref = 0, currSuff = 0, totalSum = 0;

        //calculate prefix, suffix, total sum
        for (int i = 0; i < n; i++)
        {
            currPref += a[i];
            currSuff += a[n - i - 1];
            totalSum += a[i];
            maxPref = max(maxPref, currPref);
            maxSuff = max(maxSuff, currSuff);
        }

        //if totaSum>0, then only we should add all repetitions of A, otherwise it will just get smaller with every repetition
        if (totalSum > 0)
            maxSum = maxPref + totalSum * (k - 2) + maxSuff;

        cout << maxSum << endl;
    }
}

// Max Circular Subarray Sum
long long maxCircularSum(vector<int> &arr)
{
    long long minSum = arr[0], maxSum = arr[0], minSoFar = 0, maxSoFar = 0, totalSum = 0;

    //Kadannes algo to calculate minSubarray sum and maxSubarray sum
    for (int i = 0; i < arr.size(); i++)
    {
        totalSum += arr[i];
        minSoFar += arr[i];
        maxSoFar += arr[i];
        minSum = min(minSoFar, minSum);
        maxSum = max(maxSoFar, maxSum);
        if (minSoFar > 0)
            minSoFar = 0;
        if (maxSoFar < 0)
            maxSoFar = 0;
    }

    //if all -ve , then return the maxSubarray sum
    if (totalSum == minSum)
        return maxSum;

    //else return the max of maxSubarray or the circular array
    return max(maxSum, totalSum - minSum);
}

// Subarray with given sum
/*
Use 2 pointer approach
i=0,j=0
if currSum<tar -> j++
if currSum>tar -> i++
*/
vector<int> findSubarray(vector<int> &arr, int tar)
{
    //use two pointer approach
    int i = 0, j = 0, currSum = 0;
    while (j < arr.size())
    {
        currSum += arr[j];
        if (currSum > tar)
        {
            while (i <= j && currSum > tar)
                currSum -= arr[i++];
        }
        if (currSum == tar)
            break;
        j++;
    }
    //1 based indexing, so i+1,j+1
    return {i + 1, j + 1};
}
void findSubarray()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int n, tar;
        cin >> n >> tar;

        vector<int> a(n);
        for (int i = 0; i < n; i++)
            cin >> a[i];

        vector<int> ans = findSubarray(a, tar);
        if (ans[1] > n)
        {
            cout << -1 << endl;
            continue;
        }
        cout << ans[0] << " " << ans[1] << endl;
    }
}

// 724. Find Pivot Index (Equilibrium Point)
/*
Note- The Eq point itself is not part of the left or right Sum 
use leftSum, rightSum
find total Sum of array=rsum
start iterating from begin, 
lsum=sum+a[i-1]
rsum=rsum-a[i]
*/
int pivotIndex(vector<int> &a)
{
    int n = a.size();
    if (n == 0)
        return -1;
    int lsum = 0, rsum = 0; //left = 0 , right = total
    for (int i = 0; i < n; i++)
        rsum += a[i];

    if (rsum - a[0] == 0)
        return 0;

    rsum -= a[0];

    for (int i = 1; i < n; i++)
    {
        rsum -= a[i];
        lsum += a[i - 1];
        if (lsum == rsum)
        {
            return i;
        }
    }

    return -1;
}

// Convert array into Zig-Zag fashion
/*
keep a flag to check required condition nextEle>currEle or nextEle<currEle
if condition is false then swap the two elements
*/
void zigzag(vector<int> &a)
{
    int n = a.size();
    int flag = 0; //next element should be 0->inc , 1->dec

    for (int i = 0; i < n - 1; i++)
    {
        if (flag == 0 && a[i + 1] < a[i])
        {
            swap(a[i + 1], a[i]);
        }
        else if (flag == 1 && a[i + 1] > a[i])
        {
            swap(a[i + 1], a[i]);
        }
        flag ^= 1;
    }
}

// Find Pair Given Difference
/*
Sort the array
Use 2 Pointer Approach, i=0,j=1
currDiff=arr[j]-arr[i]
while(cond..){
if(currDiff>tar) i++
if(currDiff<tar) j++
} 
*/
int diffPair(vector<int> &a, int d)
{
    int n = a.size();
    sort(a.begin(), a.end());

    int diff, i = 0, j = 1;
    while (j < n && i < n)
    {
        diff = a[j] - a[i];
        //diff will always be positive, so i<=j is true always
        if (diff > d)
        {
            i++;
        }
        else if (diff < d)
        {
            j++;
        }
        else if (i != j && diff == d)
        {
            return 1;
        }
    }

    return -1;
}

// Chocolate Distribution Problem (no submission option available on G4G)
int minDiff(vector<int> &packets, int children)
{
    sort(packets.begin(), packets.end());

    int minDiff = 1e8;
    for (int i = 0; i + children - 1 < packets.size(); i++)
    {
        minDiff = min(minDiff, packets[i + children - 1] - packets[i]);
    }

    return minDiff;
}

// Minimum Number of Platforms Required for a Railway/Bus Station
/*
sort both arrival and departure times
then do it like merging two sorted arrays
*/
int getStations(vector<int> &arrival, vector<int> &departure)
{
    int n = arrival.size();
    int currStations = 0, minStations = 0;
    sort(arrival.begin(), arrival.end());
    sort(departure.begin(), departure.end());

    int i = 0, j = 0;
    //we only need to check for arrival array for size check as all trains will arrive first,
    //so arrival array will finish first always
    while (i < n && j < n)
    {
        //check for <= as arrival and departure times can be same as well and we need a seperate stations in this case
        if (arrival[i] <= departure[j])
        {
            i++;
            currStations++;
        }
        else
        {
            j++;
            currStations--;
        }
        minStations = max(currStations, minStations);
    }

    return minStations;
}
void getStations()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int n;
        cin >> n;
        vector<int> arrival(n), departure(n);
        for (int i = 0; i < n; i++)
            cin >> arrival[i];
        for (int i = 0; i < n; i++)
            cin >> departure[i];

        cout << getStations(arrival, departure) << endl;
    }
}

// 373. Find K Pairs with Smallest Sums
/*
using min heap
the approach is similar to merging K sorted Lists using Priority Queue
*/
struct compareSum
{
    bool operator()(const vector<int> &v1, const vector<int> &v2)
    {
        return v1[0] + v1[1] > v2[0] + v2[1];
    }
};
vector<vector<int>> kSmallestPairs(vector<int> &nums1, vector<int> &nums2, int k)
{
    if (nums1.size() == 0 || nums2.size() == 0 || k == 0)
        return {};

    vector<vector<int>> ans;

    priority_queue<vector<int>, vector<vector<int>>, compareSum> pq; // {nums1[i], nums2[j], j}

    for (int i = 0; i < nums1.size() && i < k; i++)
    {
        pq.push({nums1[i], nums2[0], 0});
    }

    while (k-- > 0 && pq.size() != 0)
    {
        vector<int> pr = pq.top();
        pq.pop();

        ans.push_back({pr[0], pr[1]});

        if (pr[2] == nums2.size() - 1)
            continue;

        pq.push({pr[0], nums2[pr[2] + 1], pr[2] + 1});
    }

    return ans;
}

// 26. Remove Duplicates from Sorted Array
/*
Use two pointer approach O(n)
if(a[i]==a[j]) then only j++
if(!=) a[i+1]=a[j]; 
        i++; 
        j++;
*/
int removeDuplicates(vector<int> &nums)
{
    if (nums.size() == 0)
        return 0;
    int i = 0, j = 1;
    while (j < nums.size())
    {
        if (nums[i] != nums[j])
        {
            nums[i + 1] = nums[j];
            i++;
        }
        j++;
    }

    return i + 1;
}

// 153. Find Minimum in Rotated Sorted Array
/*
use binary search
if arr[mid]<arr[hi], then smallest element i.e. pivot lies in left region i.e. region<=mid
else it lies in right region i.e. in region>mid
*/
int findMin(vector<int> &nums)
{
    int lo = 0, hi = nums.size() - 1;
    long mid;
    while (lo < hi)
    {
        mid = (lo + hi) / 2;
        if (nums[mid] < nums[hi])
            hi = mid;
        else
            lo = mid + 1;
    }

    return nums[lo];
}

// 33. Search in Rotated Sorted Array
int search_01(vector<int> &nums, int target)
{
    int lo = 0, hi = nums.size() - 1;
    long mid;
    //find smallest element(pivot)
    while (lo < hi)
    {
        mid = (lo + hi) / 2;
        if (nums[mid] < nums[hi])
            hi = mid;
        else
            lo = mid + 1;
    }

    //find which region target lies in
    int pivot = lo;
    lo=0,hi=nums.size()-1;
    if (target<=nums[hi])
    { //right half of pivot
        lo = pivot;
    }
    else
    { //left half of pivot
        hi = pivot-1;
    }

    //normal binary search in that region
    while (lo <= hi)
    {
        mid = (lo + hi) / 2;
        if (target < nums[mid])
            hi = mid - 1;
        else if (target > nums[mid])
            lo = mid + 1;
        else
            return mid;
    }
    return -1;
}
int search_02(vector<int> &nums, int target)
{
    int lo = 0, hi = nums.size() - 1;
    long mid;
    //find smallest element(pivot)
    while (lo < hi)
    {
        mid = (lo + hi) / 2;
        if (nums[mid] < nums[hi])
            hi = mid;
        else
            lo = mid + 1;
    }

    //binary search accounting for rotation
    int rot = lo, n = nums.size();
    lo = 0;
    hi = n - 1;
    while (lo <= hi)
    {
        int mid = (lo + hi) / 2;
        int realmid = (mid + rot) % n;
        if (nums[realmid] == target)
            return realmid;
        if (nums[realmid] < target)
            lo = mid + 1;
        else
            hi = mid - 1;
    }
    return -1;
}

void solve()
{
}
int main()
{
    // g++ Arrays.cpp -o out && ./out < input.txt > output.txt
    solve();
    return 0;
}