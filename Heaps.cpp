#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <queue>

using namespace std;

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

    for (int i = 0; i < nums1.size(); i++)
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

void solve()
{
}
int main()
{
    // g++ Heaps.cpp -o out && ./out < input.txt > output.txt
    solve();
    return 0;
}