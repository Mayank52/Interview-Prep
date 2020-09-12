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
            res[j][0] = min(res[j][0], intervals[i][0]);    //start will min of both starts
            res[j][1] = max(res[j][1], intervals[i][1]);    //end will be max of both ends
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

void solve()
{
}
int main()
{
    // g++ Arrays.cpp -o out && ./out < input.txt > output.txt
    solve();
    return 0;
}