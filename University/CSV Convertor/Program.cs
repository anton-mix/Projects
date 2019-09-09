using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Web;
enum oldF {
	id = 0,
	empty = 1,
	shotDescr = 2,
	name = 3,
	descr = 4,
	path = 5,
	image = 6,
	size = 7,
	count = 8,
	prise = 9
}

enum newF {
	id = 0,
	prise = 1,
	name = 2,
	shotDescr = 3,
	discr = 4,
	path = 5,
	image = 6,
	attrGroup = 7,
	count = 8,
	deltaPrise = 9
}

namespace prob {
	class Program {
		static void Main(string[] args) {
			FileStream fstream;
			byte[] array;
			string textFromFile, newFile = "",newEntry = "",id = "",buffPrise,bufPath;
			string[] LineFile = null,newPathLine,path, oldEntry,newPath, filecsv;
			string Attr = "Артикул#Цена#Название продукта#Краткое описание#Описание#Категория#Изображения товара#Attribute Group: count#Количество комбинаций#Цена Комбинации\n";
			float prise = 0,prob = 0;
			bool prov = true;
			int a = 0;
			List<string> size = new List<string>();
			List<string> buff = new List<string>();
			List<string> LineFile2 = new List<string>();

			filecsv = Directory.GetFiles("csvs");

			for (int i = 0; i < filecsv.Length; i++) {
				fstream = File.OpenRead(filecsv[i]);
				array = new byte[fstream.Length];
				fstream.Read(array, 0, array.Length);
				fstream.Close();
				textFromFile = System.Text.Encoding.Default.GetString(array);
				LineFile = textFromFile.Split('\n');
				for (int j = 0; j < LineFile.Length; j++) {
					LineFile2.Add(LineFile[j]);
				}
			}
			fstream = File.OpenRead("NewPath.txt");
			array = new byte[fstream.Length];
			fstream.Read(array,0,array.Length);
			fstream.Close();
			textFromFile = System.Text.Encoding.Default.GetString(array);
			newPathLine = textFromFile.Split('\n');

			for (int i = 0; i < LineFile2.Count; i++) {
				oldEntry = LineFile2[i].Split('#');
				
				//Исправление ошибок
				if(oldEntry[(int)oldF.size] != "" && (oldEntry[(int)oldF.size].Split(' ')[1] == "коп" || oldEntry[(int)oldF.size].Split(' ')[1] == "р")) {
					oldEntry[(int)oldF.size] = oldEntry[(int)oldF.size].Split(' ')[0] + " шт";
				}

				if(!File.Exists("image//" + oldEntry[(int)oldF.image])) {
					oldEntry[(int)oldF.image] = "4.jpg";
					a++;
				}
				path = oldEntry[(int)oldF.path].Split('|');
				bool prov2 = true;
				for (int j =0;j< newPathLine.Length;j++) {
					newPath = newPathLine[j].Split('|');
					newPath[newPath.Length - 1] = newPath[newPath.Length-1].Remove(newPath[newPath.Length - 1].Length - 1);
					if (newPath[newPath.Length-1] == path[path.Length-1]) {
						path = newPath;
						prov2 = false;
						break;
					}
				}
				if (prov2) {
					bool prov3 = true;
					for (int k = 0; k < buff.Count; k++)
						if (oldEntry[(int)oldF.path] == buff[k]) {
							prov3 = false;
							break;
						}
					if(prov3)
						buff.Add(oldEntry[(int)oldF.path]);
				}

				//Создание деректорий
				path = oldEntry[(int)oldF.path].Split('|');
				bufPath = path[0];
				for (int j = 1;j < path.Length;j++) {
					if (!Directory.Exists(bufPath)) {
						Directory.CreateDirectory(bufPath);
					}
					if(j < path.Length-1)
						bufPath = bufPath + "\\" + path[j];
				}
				//Создания файла
				if (newFile != path[path.Length - 1]) {
					if (fstream != null)
						fstream.Close();
					newFile = path[path.Length - 1];
					fstream = new FileStream(bufPath + "\\" + path[path.Length-1] + ".csv",FileMode.OpenOrCreate);
					array = System.Text.Encoding.Default.GetBytes(Attr);
					fstream.Write(array,0,array.Length);
				}
				if (oldEntry[(int)oldF.prise].Length == 0 || oldEntry[(int)oldF.count].Length == 0 || oldEntry[(int)oldF.size].Length == 0) {
					prob++;
					continue;
				}

				for (int j = 0; j < size.Count; j++) {
					if (size[j] == oldEntry[(int)oldF.size].Split(' ')[1]) {
						prov = false;
					}
				}
				if (prov) {
					size.Add(oldEntry[(int)oldF.size].Split(' ')[1]);
				}
				prov = true;

				if (oldEntry[(int)oldF.prise].Length == 0 || oldEntry[(int)oldF.count].Length == 0 || oldEntry[(int)oldF.size].Length == 0) {
					prob++;
					continue;
				}
				/*
				// Удаление "Главной"
				oldEntry[(int)oldF.path] = path[1];
				for (int j = 2;j< path.Length;j++) {
					oldEntry[(int)oldF.path] = oldEntry[(int)oldF.path] + "|" + path[j];
				}**/
				//Создание записи
				if (id != oldEntry[(int)oldF.id]) {
					id = oldEntry[(int)oldF.id];
					buffPrise = oldEntry[(int)oldF.prise].Replace('.',',');
					buffPrise = buffPrise.Remove(buffPrise.Length - 1);
					prise = float.Parse(buffPrise);
					newEntry = oldEntry[(int)oldF.id] + "#" + oldEntry[(int)oldF.prise].Remove(oldEntry[(int)oldF.prise].Length - 1) + "#" + oldEntry[(int)oldF.name] + "#" + oldEntry[(int)oldF.shotDescr] + "#" + oldEntry[(int)oldF.descr] + "#" + oldEntry[(int)oldF.path] + "#" + oldEntry[(int)oldF.image] + "#" + oldEntry[(int)oldF.size] + "#" + oldEntry[(int)oldF.count] + "#0\n";
					/*if (oldEntry[(int)oldF.size].Split(' ')[1] == "гр" || oldEntry[(int)oldF.size].Split(' ')[1] == "кг") {
						newEntry = oldEntry[(int)oldF.id] + "#" + oldEntry[(int)oldF.prise].Remove(oldEntry[(int)oldF.prise].Length - 1) + "#" + oldEntry[(int)oldF.name] + "#" + oldEntry[(int)oldF.shotDescr] + "#" + oldEntry[(int)oldF.descr] + "#" + oldEntry[(int)oldF.path] + "#" + oldEntry[(int)oldF.image] + "##" + oldEntry[(int)oldF.size] + "###" + oldEntry[(int)oldF.count] + "#0\n";
					}
					else if (oldEntry[(int)oldF.size].Split(' ')[1] == "мл" || oldEntry[(int)oldF.size].Split(' ')[1] == "л") {
						newEntry = oldEntry[(int)oldF.id] + "#" + oldEntry[(int)oldF.prise].Remove(oldEntry[(int)oldF.prise].Length - 1) + "#" + oldEntry[(int)oldF.name] + "#" + oldEntry[(int)oldF.shotDescr] + "#" + oldEntry[(int)oldF.descr] + "#" + oldEntry[(int)oldF.path] + "#" + oldEntry[(int)oldF.image] + "###" + oldEntry[(int)oldF.size] + "##" + oldEntry[(int)oldF.count] + "#0\n";
					}
					else if (oldEntry[(int)oldF.size].Split(' ')[1] == "шт" || oldEntry[(int)oldF.size].Split(' ')[1] == "лист") {
						newEntry = oldEntry[(int)oldF.id] + "#" + oldEntry[(int)oldF.prise].Remove(oldEntry[(int)oldF.prise].Length - 1) + "#" + oldEntry[(int)oldF.name] + "#" + oldEntry[(int)oldF.shotDescr] + "#" + oldEntry[(int)oldF.descr] + "#" + oldEntry[(int)oldF.path] + "#" + oldEntry[(int)oldF.image] + "####" + oldEntry[(int)oldF.size] + "#" + oldEntry[(int)oldF.count] + "#0\n";
					}
					else {
						newEntry = oldEntry[(int)oldF.id] + "#" + oldEntry[(int)oldF.prise].Remove(oldEntry[(int)oldF.prise].Length - 1) + "#" + oldEntry[(int)oldF.name] + "#" + oldEntry[(int)oldF.shotDescr] + "#" + oldEntry[(int)oldF.descr] + "#" + oldEntry[(int)oldF.path] + "#" + oldEntry[(int)oldF.image] + "#" + oldEntry[(int)oldF.size] + "####" + oldEntry[(int)oldF.count] + "#0\n";
					}*/
				}
				else {
					newEntry = "#######" + oldEntry[(int)oldF.size] + "#" + oldEntry[(int)oldF.count] + "#" + Convert.ToString(float.Parse(oldEntry[(int)oldF.prise].Replace('.',',')) - prise) + "\n";
					/*if (oldEntry[(int)oldF.size].Split(' ')[1] == "гр" || oldEntry[(int)oldF.size].Split(' ')[1] == "кг") {
						newEntry = "########" + oldEntry[(int)oldF.size] + "###" + oldEntry[(int)oldF.count] + "#" + Convert.ToString(float.Parse(oldEntry[(int)oldF.prise].Replace('.',',')) - prise) + "\n";
					}
					else if (oldEntry[(int)oldF.size].Split(' ')[1] == "мл" || oldEntry[(int)oldF.size].Split(' ')[1] == "л") {
						newEntry = "#########" + oldEntry[(int)oldF.size] + "##" + oldEntry[(int)oldF.count] + "#" + Convert.ToString(float.Parse(oldEntry[(int)oldF.prise].Replace('.',',')) - prise) + "\n";
					}
					else if (oldEntry[(int)oldF.size].Split(' ')[1] == "шт" || oldEntry[(int)oldF.size].Split(' ')[1] == "лист") {
						newEntry = "##########" + oldEntry[(int)oldF.size] + "#" + oldEntry[(int)oldF.count] + "#" + Convert.ToString(float.Parse(oldEntry[(int)oldF.prise].Replace('.',',')) - prise) + "\n";
					}
					else {
						newEntry = "#######" + oldEntry[(int)oldF.size] + "####" + oldEntry[(int)oldF.count] + "#" + Convert.ToString(float.Parse(oldEntry[(int)oldF.prise].Replace('.',',')) - prise) + "\n";
					}*/
				}
				array = System.Text.Encoding.Default.GetBytes(newEntry);
				fstream.Write(array,0,array.Length);
			}
			fstream.Close();
		}
	}
}
