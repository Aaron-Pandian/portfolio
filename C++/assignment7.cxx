#include <cmath>

#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <memory>
using std::make_shared;
using std::shared_ptr;

#include <vector>
using std::vector;

#include <cassert>

#include <bits/stdc++.h>

#include <list>

class Node
{
private:
public:
    int data;
    Node *next;
};

class List
{
private:
    shared_ptr<Node> head{nullptr};

public:
    List(){};
    auto headnode() { return head; };
};

Node *newNode(int new_data)
{
    Node *new_node = new Node();
    new_node->data = new_data;
    new_node->next = NULL;

    return new_node;
}

void sortedInsert(Node **head_ref, Node *new_node)
{
    Node *current;
    if (*head_ref == NULL || (*head_ref)->data >= new_node->data)
    {
        new_node->next = *head_ref;
        *head_ref = new_node;
    }
    else
    {
        current = *head_ref;
        while (current->next != NULL && current->next->data < new_node->data)
        {
            current = current->next;
        }
        new_node->next = current->next;
        current->next = new_node;
    }
}

int main()
{
    Node *head = NULL;
    int input;
    cin >> input;
    Node *new_node = newNode(input);
    sortedInsert(&head, new_node);
    while (true)
    {
        if (input == 0)
        {
            break;
        }
        else
        {
            Node *temp = head;
            vector<int> values;
            while (temp != NULL)
            {
                // cout << temp->data << ", ";
                values.push_back(temp->data);
                temp = temp->next;
            }
            for (int i = 0; i < values.size(); i++)
            {
                if (i != 0)
                {
                    cout << ", ";
                }
                cout << values[i];
            }
            cout << endl;

            cin >> input;
            new_node = newNode(input);
            sortedInsert(&head, new_node);
        }
    }

    return 0;
}