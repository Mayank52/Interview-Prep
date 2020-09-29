#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <map>
#include <algorithm>

using namespace std;

// Group words with same set of characters
vector<vector<string>> groupWords(vector<string> &words)
{
    map<vector<bool>, vector<string>> mp; //char freq: words

    //make a map of (char present:words)
    for (int i = 0; i < words.size(); i++)
    {
        vector<bool> freq(26, 0);
        for (char c : words[i])
            freq[c - 'a']=true;

        mp[freq].push_back(words[i]);
    }

    vector<vector<string>> res;

    //push all words with same same character set into the result
    for (auto ele : mp)
        res.push_back(ele.second);

    for (vector<string> &word : res)
    {
        for (string str : word)
            cout << str << " ";
        cout << endl;
    }

    return res;
}


void solve()
{
    // vector<string> words{"may", "student", "students", "dog",
    //                      "studentssess", "god", "cat", "act",
    //                      "tab", "bat", "flow", "wolf", "lambs",
    //                      "amy", "yam", "balms", "looped",
    //                      "poodle"};

    // groupWords(words);
}
int main()
{
    solve();
    return 0;
}