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
            freq[c - 'a'] = true;

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

// 1347. Minimum Number of Steps to Make Two Strings Anagram
int minSteps(string s, string t)
{
    vector<int> freq(26);
    //count the  freq of characters in str1
    for (int i = 0; i < s.size(); i++)
        freq[s[i] - 'a']++;

    int count = 0;
    //check for the characters in str2, if they have same freq
    for (int i = 0; i < t.size(); i++)
    {
        if (freq[t[i] - 'a'] == 0)
            count++;
        else
            freq[t[i] - 'a']--;
    }

    return count;
}

// 890. Find and Replace Pattern
bool match(string &word, string &pattern)
{
    unordered_map<char, char> mp;

    //1.check that all characters of word are mapped to unique values
    //map all values in word to pattern
    for (int i = 0; i < word.size(); i++)
    {
        //if it present in map but mapped to some other value then return false
        if (mp.find(word[i]) != mp.end())
        {
            if (mp[word[i]] != pattern[i])
                return false;
        }
        //else add it to the map
        mp[word[i]] = pattern[i];
    }

    //2. check if all values are mapped to unique keys
    //check if more than one value has been mapped to same value i.e. more than one key has same value
    vector<bool> seen(26, false);
    for (auto ele : mp)
    {
        //if the value has already been mapped, then return false
        if (seen[ele.second - 'a'])
            return false;
        //seeing it for first time
        seen[ele.second - 'a'] = true;
    }

    return true;
}
vector<string> findAndReplacePattern(vector<string> &words, string pattern)
{
    vector<string> res;
    for (string &word : words)
    {
        if (match(word, pattern))
            res.push_back(word);
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