#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <queue>
#include <unordered_map>

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

// 154. Find Minimum in Rotated Sorted Array II (Duplicates allowed)
/*
if arr[mid]==arr[hi], then we have no choice but to do a linear search in whole array,
as we cannot decide which part pivot lies in.
So in that case we just reduce the hi by 1
*/
int findMin(vector<int> &nums)
{
    //Worst case: not rotated array, O(n)
    int lo = 0, hi = nums.size() - 1;
    long mid;
    while (lo < hi)
    {
        mid = (lo + hi) / 2;
        if (nums[mid] < nums[hi])
            hi = mid;
        else if (nums[mid] > nums[hi])
            lo = mid + 1;
        else
        {
            /*if (i-1)th is greater that means (i)th will be the pivot
            eg- 1 1 1 1 2 1 1
            here at i=5 , a[i-1]>a[i] and i is the pivot index
            */
            if (nums[hi - 1] > nums[hi])
                return nums[hi];
            hi--;
        }
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
    lo = 0, hi = nums.size() - 1;
    if (target <= nums[hi])
    { //right half of pivot
        lo = pivot;
    }
    else
    { //left half of pivot
        hi = pivot - 1;
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

// Given a sorted and rotated array, find if there is a pair with a given sum(Not available for submission)
bool pairSum(vector<int> &nums, int target)
{
    //find pivot
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

    //use 2 pointer method(meet in the middle) using mod to keep index in range
    int pivot = lo, n = nums.size();
    int i = lo, j = pivot - 1;
    //i goes pivot -> n-1 , j goes pivot-1 -> 0
    while (i != j)
    {
        if (nums[i] + nums[j] < target)
        {
            //to keep i++ in range
            i = (i + 1) % n;
        }
        else if (nums[i] + nums[j] > target)
        {
            //to keep j-- in range
            j = (n + j - 1) % n;
        }
        else
            return true;
    }

    return false;
}

// 525. Contiguous Array
/*
we do count-- for 0s, and count++ for 1s
keep a map for count:index values
two index where count is equal, will have equal number of 0s and 1s, so at each index we just check 
if this count value is present in map and update the max subarray length
*/
int findMaxLength(vector<int> &nums)
{
    int count = 0, maxlen = 0;
    unordered_map<int, int> mp; //{count:index}
    mp[0] = -1;
    for (int i = 0; i < nums.size(); i++)
    {
        if (nums[i] == 0)
            count--;
        else
            count++;
        if (mp.find(count) != mp.end())
            maxlen = max(maxlen, i - mp[count]);
        else
            mp[count] = i;
    }
    return maxlen;
}

// 121. Best Time to Buy and Sell Stock
/*
Just find the max diff in array
Keep a min so far for each element, and update the maxProfit with the max diff i.e arr[i]-minSoFar
*/
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;
    int maxp = 0, minSoFar = prices[0];
    for (int i = 0; i < prices.size(); i++)
    {
        minSoFar = min(minSoFar, prices[i]);
        maxp = max(maxp, prices[i] - minSoFar);
    }

    return maxp;
}

// 122. Best Time to Buy and Sell Stock II
/*
Approach 1-
while the price keeps increasing keep adding it to profit,
when it decreases, dont include it
*/
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;
    int totalProfit = 0;
    for (int i = 0; i < prices.size() - 1; i++)
    {
        if (prices[i + 1] > prices[i])
            totalProfit += (prices[i + 1] - prices[i]);
    }
    return totalProfit;
}
/*
Approach 2-
keep a currMin (for the current rise in price),
the moment the price drops, you sell the stock and add arr[i]-currMin into the total profit
*/
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;
    int currMin = prices[0], totalProfit = 0;
    for (int i = 0; i < prices.size() - 1; i++)
    {
        if (prices[i] > prices[i + 1])
        {
            totalProfit += prices[i] - currMin;
            currMin = prices[i + 1];
        }
    }
    totalProfit += prices[prices.size() - 1] - currMin;

    return totalProfit;
}

// 123. Best Time to Buy and Sell Stock III
/*Approach 1-
Make a prefix and suffix array
Prefix- Stock is sold on this day
Suffix- Stock is bought on this day

Then the max profit the max sum of prefix and suffix arrays
*/
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;

    int n = prices.size();
    vector<int> prefix(n), suffix(n);

    //we can maintain the max Profit so far seperatly,
    /*
    int minBuying = prices[0], maxProfitSoFar = 0;
    for (int i = 0; i < n; i++)
    {
        minBuying = min(prices[i], minBuying);
        maxProfitSoFar = max(maxProfitSoFar, prices[i] - minBuying);
        prefix[i] = maxProfitSoFar;
    }

    int maxSelling = prices[n - 1];
    maxProfitSoFar = 0;
    for (int i = n - 1; i >= 0; i--)
    {
        maxSelling = max(maxSelling, prices[i]);
        maxProfitSoFar = max(maxProfitSoFar, maxSelling - prices[i]);
        suffix[i] = maxProfitSoFar;
    }
    */

    //or we can just use the value in i-1 index, as it already has the max value upto that point
    int minBuying = prices[0];
    prefix[0] = 0;
    for (int i = 1; i < n; i++)
    {
        minBuying = min(prices[i], minBuying);
        prefix[i] = max(prefix[i - 1], prices[i] - minBuying);
    }

    int maxSelling = prices[n - 1];
    suffix[n - 1] = 0;
    for (int i = n - 2; i >= 0; i--)
    {
        maxSelling = max(maxSelling, prices[i]);
        suffix[i] = max(suffix[i + 1], maxSelling - prices[i]);
    }

    int maxp = 0;
    for (int i = 0; i < n; i++)
    {
        maxp = max(maxp, prefix[i] + suffix[i]);
    }

    return maxp;
}

//Approach 2 (faster)-
int maxProfit(vector<int> &prices)
{
    int sell1 = 0, sell2 = 0, buy1 = 1e8, buy2 = 1e8;
    for (int i = 0; i < prices.size(); i++)
    {
        buy1 = min(buy1, prices[i]);
        sell1 = max(sell1, prices[i] - buy1);
        buy2 = min(buy2, prices[i] - sell1);
        sell2 = max(sell2, prices[i] - buy2);
    }
    return sell2;
}

// 188. Best Time to Buy and Sell Stock IV (Not Complete)
int maxp; //max profit
//shares==1 -> have an extra share, so can only sell
//shares==0 -> dont have a share so can ony buy
//Recursion
//void type
void maxProfit_rec1(int shares, int profit, int K, int idx, vector<int> &prices)
{
    maxp = max(maxp, profit);
    if (K == 0 || idx == prices.size())
    {
        return;
    }

    //buy
    if (shares == 0)
    {
        maxProfit_rec1(1, profit - prices[idx], K, idx + 1, prices);
    }

    //sell
    else
    {
        maxProfit_rec1(0, profit + prices[idx], K - 1, idx + 1, prices);
    }

    //do nothing
    maxProfit_rec1(shares, profit, K, idx + 1, prices);
}

//return type
int maxProfit_rec2(int shares, int profit, int K, int idx, vector<int> &prices)
{

    if (K == 0 || idx == prices.size())
    {
        return profit;
    }

    int maxp = 0;
    //buy
    if (shares == 0)
    {
        maxp = max(maxp, maxProfit_rec2(1, profit - prices[idx], K, idx + 1, prices));
    }

    //sell
    else
    {
        maxp = max(maxp, maxProfit_rec2(0, profit + prices[idx], K - 1, idx + 1, prices));
    }

    //do nothing
    maxp = max(maxp, maxProfit_rec2(shares, profit, K, idx + 1, prices));

    return maxp;
}
int maxProfit(int k, vector<int> &prices)
{
    return maxProfit_rec2(0, 0, k, 0, prices);
}

// 714. Best Time to Buy and Sell Stock with Transaction Fee
int maxProfit(vector<int> &prices, int fee)
{
    if (prices.size() == 0)
        return 0;

    int maxAfterBuy = -prices[0]; //current cash in hand after buying
    int maxAfterSell = 0;         //current cash in hand after selling
    for (int i = 1; i < prices.size(); i++)
    {
        maxAfterBuy = max(maxAfterBuy, maxAfterSell - prices[i]);
        maxAfterSell = max(maxAfterSell, maxAfterBuy + prices[i] - fee);
    }

    return maxAfterSell;
}

// 309. Best Time to Buy and Sell Stock with Cooldown
int maxProfit(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;

    int maxAfterBuy = -prices[0]; //current cash in hand after buying
    int maxAfterSell = 0;         //current cash in hand after selling
    int prevMaxSell = 0;          //max cash after cooldown(stores the max sell just before cooldown)
    for (int i = 1; i < prices.size(); i++)
    {
        int prevMaxBuy = maxAfterBuy;
        maxAfterBuy = max(maxAfterBuy, prevMaxSell - prices[i]);
        prevMaxSell = maxAfterSell;
        maxAfterSell = max(maxAfterSell, prevMaxBuy + prices[i]);
    }

    return maxAfterSell;
}

// Maximum Difference (given that second element is greater than first element)
int maxDiff(vector<int> &prices)
{
    if (prices.size() == 0)
        return 0;
    int maxp = 0, minSoFar = prices[0];
    for (int i = 0; i < prices.size(); i++)
    {
        minSoFar = min(minSoFar, prices[i]);
        maxp = max(maxp, prices[i] - minSoFar);
    }

    //if no secondEle > firstEle then return -1
    return maxp == 0 ? -1 : maxp;
}
int maxDiff()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int n;
        cin >> n;
        vector<int> arr(n);
        for (int i = 0; i < n; i++)
            cin >> arr[i];
        cout << maxDiff(arr) << endl;
    }
}

// 56. Merge Intervals
vector<vector<int>> merge(vector<vector<int>> &intervals)
{
    if (intervals.size() == 0)
        return {};

    vector<vector<int>> res;
    //sort the array
    sort(intervals.begin(), intervals.end());
    res.push_back(intervals[0]);
    int j = 0;

    //for each interval check if it overlaps with last one, and add it accordingly
    for (int i = 1; i < intervals.size(); i++)
    {
        //if it overlaps, merge the two intervals
        if (intervals[i][0] <= res[j][1])
        {
            res[j][0] = min(res[j][0], intervals[i][0]); //start will min of both starts
            res[j][1] = max(res[j][1], intervals[i][1]); //end will be max of both ends
        }
        //if it does not overlap, just add it
        else
        {
            res.push_back(intervals[i]);
            j++;
        }
    }

    return res;
}

// 75. Sort Colors (3 Way Partition)
void sortColors(vector<int> &nums)
{
    int lt = 0, gt = nums.size() - 1, i = 0;
    while (i <= gt)
    {
        //left region: i++, lt++ (after swapping element will definitely belong to left region)
        if (nums[i] < 1)
        {
            swap(nums[lt++], nums[i++]);
        }
        //right region: only gt-- (because after swapping we could get an element that may belong to left or middle region)
        else if (nums[i] > 1)
        {
            swap(nums[gt--], nums[i]);
        }
        //middle region: i++
        else
            i++;
    }
}

// Three way partitioning
vector<int> threeWayPartition(vector<int> nums, int a, int b)
{
    int lt = 0, gt = nums.size() - 1, i = 0;
    while (i <= gt)
    {
        if (nums[i] < a)
        {
            swap(nums[lt++], nums[i++]);
        }
        else if (nums[i] > b)
        {
            swap(nums[gt--], nums[i]);
        }
        else
            i++;
    }

    return nums;
}

// 912. Sort an Array
void mergeSort(vector<int> &arr, int si, int ei)
{
    if (si == ei)
        return;

    long mid = (si + ei) / 2;
    mergeSort(arr, si, mid);
    mergeSort(arr, mid + 1, ei);

    merge(arr, si, mid, ei);
}
void merge(vector<int> &arr, int si, int mid, int ei)
{
    //make a temp res array, to store sorted list, then copy it into original array
    //size of res= ei-si+1, but it will be initialized by zero if you provide size,
    // so will have to use an index instead of push_back
    vector<int> res;
    int i = si, j = mid + 1;     //start index of both halfs
    int n = mid + 1, m = ei + 1; //end index of both halfs

    while (i < n && j < m)
    {
        if (arr[i] < arr[j])
            res.push_back(arr[i++]);
        else
            res.push_back(arr[j++]);
    }

    while (i < n)
        res.push_back(arr[i++]);
    while (j < m)
        res.push_back(arr[j++]);

    for (int k = si; k < ei + 1; k++)
        arr[k] = res[k - si];
}
vector<int> sortArray(vector<int> &nums)
{
    mergeSort(nums, 0, nums.size() - 1);
    return nums;
}

// Inversion of array
long merge_(vector<int> &arr, int si, int mid, int ei)
{
    vector<int> res;
    int i = si, j = mid + 1;
    int n = mid + 1, m = ei + 1;
    long count = 0;

    while (i < n && j < m)
    {
        if (arr[i] <= arr[j])
            res.push_back(arr[i++]);
        else
        {
            res.push_back(arr[j++]);
            //when the element in first subarray > element in second subarray, then there is an inversion
            count += n - i;
        }
    }

    while (i < n)
        res.push_back(arr[i++]);
    while (j < m)
        res.push_back(arr[j++]);

    for (int k = si; k < ei + 1; k++)
        arr[k] = res[k - si];

    return count;
}
long mergeSort_(vector<int> &arr, int si, int ei)
{
    if (si == ei)
        return 0;

    long mid = (si + ei) / 2;
    long lcount = mergeSort_(arr, si, mid);
    long rcount = mergeSort_(arr, mid + 1, ei);

    long myCount = merge_(arr, si, mid, ei);

    return lcount + rcount + myCount;
}
long countInversions(vector<int> &arr)
{
    return mergeSort_(arr, 0, arr.size() - 1);
}

// 775. Global and Local Inversions
/*
Local inverions will occur between consecutive elements only, 
so if there is a single inversion not between consecutive elements, 
then there will be more global than local inversions.

All local inversions are global inversions.
If the number of global inversions is equal to the number of local inversions,
it means that all global inversions in permutations are local inversions.
It also means that we can not find A[i] > A[j] with i+2<=j.
meaning for global==local to be true, inversions can only be between i,i+1
*/
bool isIdealPermutation(vector<int> &A)
{
    if (A.size() < 2)
        return true;
    int currMax = 0;
    for (int i = 0; i < A.size() - 2; i++)
    {
        currMax = max(currMax, A[i]);
        //if at any point current max element is greater than the i+2 position element, then global>local
        if (currMax > A[i + 2])
            return false;
    }
    return true;
}

// 560. Subarray Sum Equals K
int subarraySum(vector<int> &nums, int k)
{
    unordered_map<int, int> mp; //prefixSum : count (how many times it has occured till now)
    int currSum = 0, count = 0;

    for (int i = 0; i < nums.size(); i++)
    {
        currSum += nums[i];
        if (currSum == k)
            count++;
        if (mp.find(currSum - k) != mp.end())
            count += mp[currSum - k];
        mp[currSum]++;
    }

    return count;
}

// 152. Maximum Product Subarray
int maxProduct(vector<int> &nums)
{
    if (nums.size() == 0)
        return 0;

    int maxpos = nums[0], minneg = nums[0], omax = nums[0]; //omax=overall max

    for (int i = 1; i < nums.size(); i++)
    {
        if (nums[i] < 0)
            swap(minneg, maxpos);

        maxpos = max(nums[i], maxpos * nums[i]);
        minneg = min(nums[i], minneg * nums[i]);

        omax = max(omax, maxpos);
    }

    return omax;
}

// Minimize the heights
/*
Just remember it (not sure about logic)
*/
int getMinDiff(int arr[], int n, int k)
{
    //sort the heights
    sort(arr, arr + n);

    //find min difference in before changing with +K or -K
    int minDiff = arr[n - 1] - arr[0];

    int minEle = arr[0] + k, maxEle = arr[n - 1] - k;
    if (maxEle < minEle)
        swap(maxEle, minEle);

    //for the whole array, find the min and max you can get
    for (int i = 1; i < n - 1; i++)
    {
        int currLargest = arr[i] + k;
        int currSmallest = arr[i] - k;

        //if curr Ele does not become both currMax and currMin after +k and -k then skip it
        if (currLargest < maxEle || currSmallest > minEle)
            continue;

        //choose whether it gives minDiff after +k or -k
        if (maxEle - currSmallest <= currLargest - minEle)
            minEle = currSmallest;
        else
            maxEle = currLargest;
    }

    //return min of the original mindiff and after +K and -K
    return min(minDiff, maxEle - minEle);
}

// Minimum swaps and K together
/*
Using Sliding window
Find count of all elements which are less than or equals to ‘k’. Let’s say the count is ‘cnt’
Using two pointer technique for window of length ‘cnt’, each time keep track of how many elements in this range are greater than ‘k’.
Let’s say the total count is ‘bad’.
Repeat step 2, for every window of length ‘cnt’ and take minimum of count ‘bad’ among them. This will be the final answer.
*/
long minSwaps(vector<int> &arr, int k)
{
    int n = arr.size();
    long count = 0, currCount = 0, minCount;

    //find total count
    for (int i = 0; i < n; i++)
        if (arr[i] <= k)
            count++;

    //find count for first window
    for (int i = 0; i < count; i++)
    {
        if (arr[i] > k)
            currCount++;
    }

    //after that, just for each window, just check the element that was removed from it (i-count), and added to it(current i)
    minCount = currCount;
    for (int i = count; i < n; i++)
    {
        if (arr[i - count] > k)
            currCount--;
        if (arr[i] > k)
            currCount++;

        minCount = min(minCount, currCount);
    }

    return minCount;
}

// Maximum sum of i*arr[i] among all rotations of a given array
/*
We can take rotations in either direction, I am taking anti clockwise(towards left)
Brute Force- O(n^2)
Calculate all sums for all rotations of array

Efficient-
Calculate sum of all arr[i]*i for original array, total Sum of all elements
after each rotation the change in sum will be due to two things:;
    the 0th element goes to n-1, arr[i]*0 = 0 becomes arr[i]*(n-1)>0 , so we do rotatedSum+=arr[i]*(n-1)
    every other element goes from i to i-1, so we subtract the total Sum from rotatedSum except for the arr[i] element 
so rotatedSum=previousRotatedSum + arr[i]*(n-1) - (totalSum-arr[i])
*/
int max_sum(int arr[], int n)
{
    long long totalSum = 0, rsum = 0, maxSum = -1e8;
    //calculate initial totalSum, rotatedSum
    for (int i = 0; i < n; i++)
    {
        totalSum += arr[i];
        rsum += arr[i] * i;
    }

    //for each rotation, find the new rotatedSum, by using the formula above
    for (int i = 0; i < n; i++)
    {
        rsum += arr[i] * (n - 1) - (totalSum - arr[i]);
        maxSum = max(maxSum, rsum);
    }

    return maxSum;
}

// 41. First Missing Positive
/*
Approach-
Put every arr[i]>0 at i+1 index
then iterate through array, and wherever arr[i]!=i+1, i+1 is the missing number
if none is missing then n+1 is the missng number
*/
int firstMissingPositive(vector<int> &nums)
{
    int n = nums.size();

    for (int i = 0; i < n; i++)
    {
        while (nums[i] > 0 && nums[i] < n + 1 && nums[nums[i] - 1] != nums[i])
            swap(nums[i], nums[nums[i] - 1]);
    }

    for (int i = 0; i < n; i++)
    {
        if (nums[i] != i + 1)
            return i + 1;
    }

    return n + 1;
}

// 283. Move Zeroes (Move all zeroes to end of array)
//Approach 1-
void moveZeroes_01(vector<int> &nums)
{
    int idx = 0;
    for (int i = 0; i < n; i++)
    {
        if (nums[i] != 0)
            nums[idx++] = nums[i];
    }
    for (int i = idx; i < nums.size(); i++)
        nums[i] = 0;
}
//Approach 2- (Optimal)
void moveZeroes_02(vector<int> &nums)
{
    /*
    instead of going through array again and changing to 0,
    we just do it during first iteration, by swapping the non zero element with last position of non zero element
    */
    int idx = 0;
    for (int i = 0; i < nums.size(); i++)
    {
        if (nums[i] != 0)
            swap(nums[idx++], nums[i]);
    }
}

// Largest Sum Subarray of Size at least K
long long largestSum(vector<int> &arr, int k)
{
    int n = arr.size();
    vector<long long> prefixMaxSum(n);

    //make prefix array of all max contiguos sum till that index
    long long maxSoFar = arr[0], currMax = 0;
    for (int i = 0; i < n; i++)
    {
        currMax += arr[i];
        maxSoFar = max(maxSoFar, currMax);
        if (currMax < 0)
            currMax = 0;
        prefixMaxSum[i] = currMax;
    }

    //for each window (i to j),find window sum, and update overall max with max(window sum, windowSum + maxPrefixSum) for the last index(i-1)
    long long currSum = 0, maxSum;
    for (int i = 0; i < k; i++)
        currSum += arr[i];
    maxSum = currSum;
    for (int i = k; i < n; i++)
    {
        //currSum = window sum
        currSum += arr[i] - arr[i - k];
        //prefixMaxSum[i-k] = max Subarray sum for the last index just before the current window
        maxSum = max(maxSum, max(currSum, currSum + prefixMaxSum[i - k]));
    }

    return maxSum;
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