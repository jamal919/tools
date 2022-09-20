// #include <config.h>
// 
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// 
// #include "tools-structures.h"
// 
// #include "list.h"
// 
// List* list_nth (List *list,int  n)
// {
//   while ((n-- > 0) && list)
//     list = list->next;
//   return list;
// }
// 
// void * list_nth_data (List*list,int n)
// {
//   while ((n-- > 0) && list)
//     list = list->next;
//   return list ? list->data : NULL;
// }
// 
// List* list_last (List *list)
// {
//   if (list) {
//       while (list->next)
//         list = list->next;
//     }
//   return list;
// }
// 
// List* list_first ( List *list)
// {
//   if (list) {
//     while (list->prev)
//       list = list->prev;
//     }
//   return list;
// }
// 
// int list_length (List *list)
// {
//   int length;
// 
//   length = 0;
//   while (list) {
//     length++;
//     list = list->next;
//     }
//   return length;
// }
// 
// List* list_append (List*list,void *data)
// {
//   List *new_list=NULL;
//   List *last=NULL;
// 
//   exitIfNull(
//     new_list=(List *)calloc(1,sizeof(List))
//     );
//   new_list->data = data;
// 
//   if (list) {
//     last = list_last (list);
//     last->next = new_list;
//     new_list->prev = last;
//     return list;
//     }
//   else
//     return new_list;
// }
// 
// List* list_prepend (List*list,void *data)
// {
//   List *new_list=NULL;
// 
//   exitIfNull(
//     new_list=(List *)calloc(1,sizeof(List))
//     );
//   new_list->data = data;
//   if (list) {
//     if (list->prev) {
//       list->prev->next = new_list;
//       new_list->prev = list->prev;
//       }
//     list->prev = new_list;
//     new_list->next = list;
//     }
// 
//   return new_list;
// }
// 
// 
// 
// List* list_insert (List*list,void *data,int position)
// {
//   List *new_list=NULL;
//   List *tmp_list=NULL;
//   
//   if (position < 0)
//     return list_append (list, data);
//   else if (position == 0)
//     return list_prepend (list, data);
//   
//   tmp_list = list_nth (list, position);
//   if (!tmp_list)
//     return list_append (list, data);
//   
//   exitIfNull(
//     new_list=(List *)calloc(1,sizeof(List))
//     );
//   new_list->data = data;
//   
//   if (tmp_list->prev) {
//       tmp_list->prev->next = new_list;
//       new_list->prev = tmp_list->prev;
//     }
//   new_list->next = tmp_list;
//   tmp_list->prev = new_list;
//   
//   if (tmp_list == list)
//     return new_list;
//   else
//     return list;
// }
// 
// List * list_free(List *list)
// {
//   List *last=NULL;
//   
//   last = list_last (list);
//   while(last!=NULL)
//     {
//       last=last->prev;
//       if(last!=NULL) {
//           free(last->next);
//           last->next=NULL;
//         }
//     }
//   free(list);
//   return(last);
// 
// }
