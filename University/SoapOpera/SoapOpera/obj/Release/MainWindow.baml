   M S B A M L   `   `   `  ���� D  @SoapOpera, Version=1.0.0.0, Culture=neutral, PublicKeyToken=null%clr-namespace:SoapOpera	SoapOpera      SoapOpera.MainWindow   R NWindowsBase, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35W SPresentationCore, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35\ XPresentationFramework, Version=4.0.0.0, Culture=neutral, PublicKeyToken=31bf3856ad364e35D 9http://schemas.microsoft.com/winfx/2006/xaml/presentation    R NSystem.Xaml, Version=4.0.0.0, Culture=neutral, PublicKeyToken=b77a5c561934e089:x,http://schemas.microsoft.com/winfx/2006/xaml     @d2http://schemas.microsoft.com/expression/blend/2008     Jmc;http://schemas.openxmlformats.org/markup-compatibility/2006     )localclr-namespace:SoapOpera       � Title$  	SoapOpera��$	��226��$	��600��.��� .+��� -    ��Name$ 
tabControl��.@��� $,��&Сформировать xlsx файл��.��� ���	#FFE5E5E5.+��� -    ��Name$ lPathFolder��$��Имя файла: ��$
��Left=�$	��Top=�$��	10,10,0,0q�$��75���� -    ��Name$ tbPahtFolder��$��23�� �� TextWrapping$
 Wrap=�$�� ��$	��Top=�$��90,10,115,0q�{�.  �� -    ��Name$ bPathFolder��$��Выбрать папку��$	��Top=�$��	0,10,10,0q�$��23��$��Right=�$	��100���� -   $ lAlbum��$��Альбомы��$
��Left=�$	��Top=�$��	10,40,0,0q�$��75���� -    ��Name$ cbAlbum��$	��Top=�$��
90,40,10,0q�$��23���� -   $ 	lDataFrom��$��Дата: c��$
��Left=�$	��Top=�$��	10,71,0,0q��� -   $ 
tbDataFrom��$
��Left=�$��23��$
 Wrap=�$�� ��$	��Top=�$	��220��$��	62,74,0,0q��� -	   $ lDataTo��$
��по��$
��Left=�$	��Top=�$��
287,71,0,0q��� -
   $ tbDataTo��$��23��$
 Wrap=�$�� ��$	��Top=�$��316,74,10,0q��� -   $ 	lProgress��$&�� Текущий прогресс:��$	��Top=�$��
10,100,0,0q�$
��Left=�$	��115��� -    �Name$ 
pbProgress��$��23��$	��Top=�$��130,103,137,0q��� -   $ 
bCreatFile��$'��!Сформировать файл��$��Right=�$	��Top=�$	��122��$��
0,104,10,0q��� -   $ lErrors��$�� ��$��10,131,10,0q��� $)��#Отправка сообщений��.��� ���	#FFE5E5E5�� $��Настройки��.��� ���	#FFE5E5E5.+��� -   $ lOptionCreate��$B��<Настройки формирования xlsx файла:��$
��Left=�$	��Top=�$��	10,10,0,0q��� -    ��Name$ cbExeptGroupId��$O��IИсключать комментарии группы из поиска.��$
��Left=�$	��Top=�$��	10,41,0,0q��� -   $ cbExeptEmptyPhoto��$U��OНе показывать фотографии без комментариев.��$
��Left=�$	��Top=�$��	10,62,0,0q�