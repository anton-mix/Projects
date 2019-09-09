using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using System.Windows.Forms;
using System.IO;
using System.Windows.Threading;
using System.Drawing;

using Microsoft.Office.Interop.Excel;
using System.Xml;


using VkNet;
using VkNet.Enums.Filters;
using VkNet.Infrastructure;
using VkNet.Model;
using VkNet.Model.RequestParams;
using System.Runtime.InteropServices;

namespace SoapOpera {
	static class MyFunctions {
		[DllImport("kernel32.dll", SetLastError = true)]
		[return: MarshalAs(UnmanagedType.Bool)]
		public static extern bool TerminateProcess(IntPtr hProcess, uint uExitCode);


		public static List<People> GetPeopleFromAlbum(int SelectedAlbum, VkNet.Utils.VkCollection<PhotoAlbum> AllAlbums, VkApi api, DateTime TimeFrom, DateTime TimeTo) {
			ulong offset = 0, count = 500;
			List<VkNet.Model.Attachments.Photo> AllPhoto = new List<VkNet.Model.Attachments.Photo>();

			try {
				while (true) {
					var photos = api.Photo.Get(new PhotoGetParams { OwnerId = AllAlbums[SelectedAlbum].OwnerId, Extended = true, Count = count, Offset = offset, AlbumId = VkNet.Enums.SafetyEnums.PhotoAlbumType.Id(AllAlbums[SelectedAlbum].Id) });
					for (int i = 0; i < photos.Count; i++) {
						AllPhoto.Add(photos[i]);
					}
					if (AllPhoto.Count == AllAlbums[SelectedAlbum].Size)
						break;
					else if (AllPhoto.Count > AllAlbums[SelectedAlbum].Size) {
						throw new OutOfBorderPhotoException();
					}
					else {
						offset = offset + count;
					}
				}
			}
			catch {
				throw new ErrorFindPhotoException(AllAlbums[SelectedAlbum].Title);
			}


			List<List<VkNet.Model.Comment>> AllComents = new List<List<VkNet.Model.Comment>>();
			for (int i = 0; i < AllPhoto.Count; i++) {
				AllComents.Add(new List<VkNet.Model.Comment>());
				if (AllPhoto[i].Comments.Count != 0) {
					count = 100;
					offset = 0;
					bool stop = true;
					try {
						while (true) {
							var Comments = api.Photo.GetComments(new PhotoGetCommentsParams { Sort = VkNet.Enums.SafetyEnums.CommentsSort.Desc, Count = count, Offset = offset, OwnerId = AllAlbums[SelectedAlbum].OwnerId, PhotoId = (ulong)AllPhoto[i].Id });
							for (int j = 0; j < Comments.Count; j++) {
								if (Comments[j].Date >= TimeFrom && Comments[j].Date <= TimeTo) {
									AllComents[i].Add(Comments[j]);
								}
								else {
									stop = false;
									break;
								}
							}
							if (AllComents[i].Count == AllPhoto[i].Comments.Count || !stop)
								break;
							else if (AllComents[i].Count > AllPhoto[i].Comments.Count) {
								throw new OutOfBorderCommentsException();
							}
							else {
								offset = offset + count;
							}
						}
					}
					catch {
						throw new ErrorFindCommentsException(i, AllAlbums[SelectedAlbum].Title);
					}
				}
			}
			List<People> peoples = new List<People>();
			for (int i = 0; i < AllComents.Count; i++) {
				for (int j = 0; j < AllComents[i].Count; j++) {
					People buff;
					if ((buff = peoples.Find(x => x.Id == AllComents[i][j].FromId)) != null) {
						buff.CommentsText[i].Add(AllComents[i][j].Text);
						buff.CommentsData[i].Add(AllComents[i][j].Date.Value);
					}
					else {
						peoples.Add(new People(AllComents.Count, "", AllComents[i][j].FromId));
						peoples[peoples.Count - 1].CommentsText[i].Add(AllComents[i][j].Text);
						peoples[peoples.Count - 1].CommentsData[i].Add(AllComents[i][j].Date.Value);
					}
				}
			}
			return peoples;
		}

		static public void WriteExcelFile(List<List<People>> peoples, List<string> NameAlbum, bool ExeptEmptyPhoto, List<List<string>> AllNamePhoto, int Format,  String NameFile) {
			Microsoft.Office.Interop.Excel.Application excel = new Microsoft.Office.Interop.Excel.Application();
			excel.Visible = false;
			excel.DisplayAlerts = false;
			Workbook worKbooK = excel.Workbooks.Add(Type.Missing);
			Worksheet worKsheeT = (Worksheet)worKbooK.ActiveSheet;
			worKsheeT.Name = "SoapOpera";
			int line = 1;
			try {
				for (int z = 0; z < peoples.Count; z++) {
					if (peoples.Count > 1 && peoples[z].Count == 0)
						continue;
					if (peoples.Count == 1 && peoples[z].Count == 0)
						throw new NotCommentException();

					bool[] buffPhoto = new bool[peoples[z][0].CommentsText.Count];
					int curentLenght = 0;
					if (ExeptEmptyPhoto) {
						for (int j = 0; j < peoples[z].Count; j++)
							for (int k = 0; k < peoples[z][0].CommentsText.Count; k++)
								if (peoples[z][j].CommentsText[k].Count > 0)
									buffPhoto[k] = true;
					}
					else {
						for (int k = 0; k < peoples[z][0].CommentsText.Count; k++)
							buffPhoto[k] = true;
					}

					for (int j = 0; j < peoples[z][0].CommentsText.Count; j++)
						if (buffPhoto[j])
							curentLenght++;

					worKsheeT.Range[worKsheeT.Cells[line, 1], worKsheeT.Cells[line + 2, 5 + curentLenght * 3]].NumberFormat = "@";
					worKsheeT.Cells[line, 1] = NameAlbum[z];

					line++;

					worKsheeT.Cells[line, 1] = "№";
					worKsheeT.Cells[line, 2] = "Участник";
					worKsheeT.Cells[line, 3] = "id";
					int curcounter = 0;
					for (int i = 0; i < peoples[z][0].CommentsText.Count; i++) {
						if (buffPhoto[i]) {
							worKsheeT.Cells[line, 4 + curcounter * 3] = AllNamePhoto[z][i].Replace('\n', ' ');
							curcounter++;
						}
					}
					worKsheeT.Cells[line, 4 + curentLenght * 3] = "Сумма к оплате";
					worKsheeT.Cells[line, 5 + curentLenght * 3] = "Предоставлено к оплате";
					line++;

					for (int i = 4; i < 3 + curentLenght * 3; i = i + 3) {
						worKsheeT.Cells[line, i] = "Комментарии";
						worKsheeT.Cells[line, i + 1] = "Кол-во";
						worKsheeT.Cells[line, i + 2] = "Цена";
					}

					line++;
					int lineStartPeople = line;
					for (int i = 0; i < peoples[z].Count; i++) {
						int offset = 0;
						string[] buff = new string[peoples[z][i].CommentsText.Count];
						for (int j = 0; j < peoples[z][i].CommentsText.Count; j++)
							buff[j] = "";
						while (true) {
							worKsheeT.Range[worKsheeT.Cells[line, 1], worKsheeT.Cells[line, 5 + curentLenght * 3]].NumberFormat = "@";
							worKsheeT.Cells[line, 1] = (i + 1);
							worKsheeT.Cells[line, 2] = peoples[z][i].Name;
							worKsheeT.Cells[line, 3] = peoples[z][i].Id;
							worKsheeT.Cells[line, 3].Interior.Color = XlRgbColor.rgbGray;
							bool prov = false;
							curcounter = 4;
							for (int j = 4; j < 4 + peoples[z][i].CommentsText.Count * 3; j = j + 3) {
								if (buffPhoto[(j - 4) / 3]) {
									worKsheeT.Range[worKsheeT.Cells[line, curcounter + 1], worKsheeT.Cells[line, curcounter + 2]].NumberFormat = "0";
									worKsheeT.Cells[line, curcounter + 2].Interior.Color = XlRgbColor.rgbYellow;
									if (peoples[z][i].CommentsText[(j - 4) / 3].Count > offset) {
										if (buff[(j - 4) / 3].Length == 0)
											buff[(j - 4) / 3] = peoples[z][i].CommentsText[(j - 4) / 3][offset];
										else
											buff[(j - 4) / 3] = buff[(j - 4) / 3] + " // " + peoples[z][i].CommentsText[(j - 4) / 3][offset];
										worKsheeT.Cells[line, curcounter] = buff[(j - 4) / 3];
										worKsheeT.Cells[line, curcounter + 1] = worKsheeT.Cells[line, curcounter + 2] = "";
										prov = true;
									}
									curcounter = curcounter + 3;
								}
							}
							if (!prov) {
								worKsheeT.Range[worKsheeT.Cells[line, curcounter], worKsheeT.Cells[line, curcounter + 1]].NumberFormat = "0";
								string buff2 = "= " + convertTo(5) + line + " * " + convertTo(6) + line;
								for (int j = 1; j < curentLenght; j++) {
									buff2 = buff2 + " + " + convertTo(5 + j * 3) + line + " * " + convertTo(6 + j * 3) + line;
								}
								worKsheeT.Cells[line, curcounter] = buff2;
								line++;
								break;
							}
							else {
								offset++;
							}
						}
					}

					worKsheeT.Range[worKsheeT.Cells[line, 1], worKsheeT.Cells[line, 5 + curentLenght * 3]].Interior.Color = Microsoft.Office.Interop.Excel.XlRgbColor.rgbWhite;

					Microsoft.Office.Interop.Excel.Borders border = worKsheeT.Range[worKsheeT.Cells[lineStartPeople, 1], worKsheeT.Cells[line - 1, 5 + curentLenght * 3]].Borders;
					border.LineStyle = Microsoft.Office.Interop.Excel.XlLineStyle.xlContinuous;
					border.Weight = 2d;


					Microsoft.Office.Interop.Excel.Range usedrange = worKsheeT.Range[worKsheeT.Cells[lineStartPeople - 3, 1], worKsheeT.Cells[lineStartPeople - 1, 5 + curentLenght * 3]];
					usedrange.Rows.AutoFit();
					usedrange.Columns.AutoFit();

					usedrange = worKsheeT.Range[worKsheeT.Cells[lineStartPeople - 3, 2], worKsheeT.Cells[line - 1, 2]];
					usedrange.Rows.AutoFit();
					usedrange.Columns.AutoFit();
					line++;
				}

				if (Format <= 0)
					worKbooK.SaveAs(NameFile, XlFileFormat.xlExcel3);
				else
					worKbooK.SaveAs(NameFile);

				worKbooK.Close();
			}
			catch (Exception e) {
				throw e;
			}
			finally {
				IntPtr handle = (IntPtr)excel.Hwnd;
				excel.Quit();
				TerminateProcess(handle, 1);
			}
		}


		static public void WriteExcelFileWithCollum(List<List<People>> peoples, List<string> NameAlbum, bool ExeptEmptyPhoto, List<List<string>> AllNamePhoto, int Format, String NameFile) {
			Microsoft.Office.Interop.Excel.Application excel = new Microsoft.Office.Interop.Excel.Application();
			excel.Visible = false;
			excel.DisplayAlerts = false;
			Workbook worKbooK = excel.Workbooks.Add(Type.Missing);
			Worksheet worKsheeT = (Worksheet)worKbooK.ActiveSheet;
			worKsheeT.Name = "SoapOpera";
			int line = 1;
			try {
				for (int z = 0; z < peoples.Count; z++) {
					if (peoples.Count > 1 && peoples[z].Count == 0)
						continue;
					if (peoples.Count == 1 && peoples[z].Count == 0)
						throw new NotCommentException();

					bool[] buffPhoto = new bool[peoples[z][0].CommentsText.Count];
					int curentLenght = 0;
					if (ExeptEmptyPhoto) {
						for (int j = 0; j < peoples[z].Count; j++)
							for (int k = 0; k < peoples[z][0].CommentsText.Count; k++)
								if (peoples[z][j].CommentsText[k].Count > 0)
									buffPhoto[k] = true;
					}
					else {
						for (int k = 0; k < peoples[z][0].CommentsText.Count; k++)
							buffPhoto[k] = true;
					}

					for (int j = 0; j < peoples[z][0].CommentsText.Count; j++)
						if (buffPhoto[j])
							curentLenght++;

					worKsheeT.Range[worKsheeT.Cells[line, 1], worKsheeT.Cells[line + 2, 5 + curentLenght * 4]].NumberFormat = "@";
					worKsheeT.Cells[line, 1] = NameAlbum[z];

					line++;

					worKsheeT.Cells[line, 1] = "№";
					worKsheeT.Cells[line, 2] = "Участник";
					worKsheeT.Cells[line, 3] = "id";
					int curcounter = 0;
					for (int i = 0; i < peoples[z][0].CommentsText.Count; i++) {
						if (buffPhoto[i]) {
							worKsheeT.Cells[line, 4 + curcounter * 4] = AllNamePhoto[z][i].Replace('\n', ' ');
							curcounter++;
						}
					}
					worKsheeT.Cells[line, 4 + curentLenght * 4] = "Сумма к оплате";
					worKsheeT.Cells[line, 5 + curentLenght * 4] = "Предоставлено к оплате";
					line++;

					for (int i = 4; i < 3 + curentLenght * 4; i = i + 4) {
						worKsheeT.Cells[line, i] = "Комментарии";
						worKsheeT.Cells[line, i + 1] = "Размерность";
						worKsheeT.Cells[line, i + 2] = "Кол-во";
						worKsheeT.Cells[line, i + 3] = "Цена";
					}

					line++;
					int lineStartPeople = line;
					for (int i = 0; i < peoples[z].Count; i++) {
						int offset = 0;
						string[] buff = new string[peoples[z][i].CommentsText.Count];
						for (int j = 0; j < peoples[z][i].CommentsText.Count; j++)
							buff[j] = "";
						while (true) {
							worKsheeT.Range[worKsheeT.Cells[line, 1], worKsheeT.Cells[line, 5 + curentLenght * 4]].NumberFormat = "@";
							worKsheeT.Cells[line, 1] = (i + 1);
							worKsheeT.Cells[line, 2] = peoples[z][i].Name;
							worKsheeT.Cells[line, 3] = peoples[z][i].Id;
							worKsheeT.Cells[line, 3].Interior.Color = XlRgbColor.rgbGray;
							bool prov = false;
							curcounter = 4;
							for (int j = 4; j < 4 + peoples[z][i].CommentsText.Count * 3; j = j + 3) {
								if (buffPhoto[(j - 4) / 3]) {

									worKsheeT.Range[worKsheeT.Cells[line, curcounter + 2], worKsheeT.Cells[line, curcounter + 3]].NumberFormat = "0";
									worKsheeT.Cells[line, curcounter + 3].Interior.Color = XlRgbColor.rgbYellow;
									if (peoples[z][i].CommentsText[(j - 4) / 3].Count > offset) {
										if (buff[(j - 4) / 3].Length == 0)
											buff[(j - 4) / 3] = peoples[z][i].CommentsText[(j - 4) / 3][offset];
										else
											buff[(j - 4) / 3] = buff[(j - 4) / 3] + " // " + peoples[z][i].CommentsText[(j - 4) / 3][offset];
										worKsheeT.Cells[line, curcounter] = buff[(j - 4) / 3];
										worKsheeT.Cells[line, curcounter + 1] = worKsheeT.Cells[line, curcounter + 2] = worKsheeT.Cells[line, curcounter + 3] = "";
										prov = true;
									}
									curcounter = curcounter + 4;
								}
							}
							if (!prov) {
								worKsheeT.Range[worKsheeT.Cells[line, curcounter], worKsheeT.Cells[line, curcounter + 1]].NumberFormat = "0";
								string buff2 = "= " + convertTo(6) + line + " * " + convertTo(7) + line;
								for (int j = 1; j < curentLenght; j++) {
									buff2 = buff2 + " + " + convertTo(6 + j * 4) + line + " * " + convertTo(7 + j * 4) + line;
								}
								worKsheeT.Cells[line, curcounter] = buff2;
								line++;
								break;
							}
							else {
								offset++;
							}
						}
					}

					worKsheeT.Range[worKsheeT.Cells[line, 1], worKsheeT.Cells[line, 5 + curentLenght * 4]].Interior.Color = Microsoft.Office.Interop.Excel.XlRgbColor.rgbWhite;

					Microsoft.Office.Interop.Excel.Borders border = worKsheeT.Range[worKsheeT.Cells[lineStartPeople, 1], worKsheeT.Cells[line - 1, 5 + curentLenght * 4]].Borders;
					border.LineStyle = Microsoft.Office.Interop.Excel.XlLineStyle.xlContinuous;
					border.Weight = 2d;


					Microsoft.Office.Interop.Excel.Range usedrange = worKsheeT.Range[worKsheeT.Cells[lineStartPeople - 3, 1], worKsheeT.Cells[lineStartPeople - 1, 5 + curentLenght * 4]];
					usedrange.Rows.AutoFit();
					usedrange.Columns.AutoFit();

					usedrange = worKsheeT.Range[worKsheeT.Cells[lineStartPeople - 3, 2], worKsheeT.Cells[line - 1, 2]];
					usedrange.Rows.AutoFit();
					usedrange.Columns.AutoFit();
					line++;
				}

				if (Format <= 0)
					worKbooK.SaveAs(NameFile, XlFileFormat.xlExcel3);
				else
					worKbooK.SaveAs(NameFile);

				worKbooK.Close();
			}
			catch (Exception e) {
				throw e;
			}
			finally {
				IntPtr handle = (IntPtr)excel.Hwnd;
				excel.Quit();
				TerminateProcess(handle, 1);
			}
		}

		static string convertTo(int number) {
			int div = number;
			string str = "";
			int mod = 0;
			while (div > 0) {
				mod = (div - 1) % 26;
				str = Convert.ToChar(65 + mod).ToString() + str;
				div = (int)((div - mod) / 26);
			}
			return str;
		}

		public static SendAlbum GetSendPeoples(string FileName) {
			SendAlbum album = new SendAlbum();
			int line = 1;
			int collCount = 4;
			Microsoft.Office.Interop.Excel.Application excel = new Microsoft.Office.Interop.Excel.Application();
			Workbook worKbooK = null;
			excel.Visible = false;
			excel.DisplayAlerts = false;


			try {
				worKbooK = excel.Workbooks.Open(FileName);
				Worksheet worKsheeT = (Worksheet)worKbooK.ActiveSheet;


				line = 1;
				album.Name = (string)(worKsheeT.Cells[line, 1] as Microsoft.Office.Interop.Excel.Range).Value;
				album.Name = (album.Name.Split('.').Length == 1) ? album.Name : album.Name.Split('.')[1];
				line = 2;
				while (true) {
					if (Convert.ToString((worKsheeT.Cells[line, collCount] as Microsoft.Office.Interop.Excel.Range).Value) == "Сумма к оплате")
						break;
					else if (Convert.ToString((worKsheeT.Cells[line, collCount] as Microsoft.Office.Interop.Excel.Range).Value) != null && Convert.ToString((worKsheeT.Cells[line, collCount] as Microsoft.Office.Interop.Excel.Range).Value) != "") 
						album.Photos.Add(Convert.ToString((worKsheeT.Cells[line, collCount] as Microsoft.Office.Interop.Excel.Range).Value));
					collCount ++;
				}
				line = 4;

				while (true) {
					if (Convert.ToString((worKsheeT.Cells[line, 1] as Microsoft.Office.Interop.Excel.Range).Value) == "" || Convert.ToString((worKsheeT.Cells[line, 1] as Microsoft.Office.Interop.Excel.Range).Value) == null) {
						break;
					}

					album.Peoples.Add(new SendPeople());
					album.Peoples[album.Peoples.Count - 1].Name = (string)(worKsheeT.Cells[line, 2] as Microsoft.Office.Interop.Excel.Range).Value;
					album.Peoples[album.Peoples.Count - 1].Id = Convert.ToInt64((worKsheeT.Cells[line, 3] as Microsoft.Office.Interop.Excel.Range).Value);
					collCount = 4;
					SendProduct buffItem = new SendProduct();
					bool breakwhile = true;
					while (breakwhile) {
						if (Convert.ToString((worKsheeT.Cells[3, collCount] as Microsoft.Office.Interop.Excel.Range).Value) == "" || Convert.ToString((worKsheeT.Cells[3, collCount] as Microsoft.Office.Interop.Excel.Range).Value) == null) {
							if (Convert.ToString((worKsheeT.Cells[2, collCount] as Microsoft.Office.Interop.Excel.Range).Value) == "Сумма к оплате") {
								album.Peoples[album.Peoples.Count - 1].AllCost = Convert.ToDouble((worKsheeT.Cells[line, collCount] as Microsoft.Office.Interop.Excel.Range).Value);
							}
							break;
						}

						switch ((string)(worKsheeT.Cells[3, collCount] as Microsoft.Office.Interop.Excel.Range).Value) {
							case "Комментарии":
								buffItem.Name = Convert.ToString((worKsheeT.Cells[2, collCount] as Microsoft.Office.Interop.Excel.Range).Value);
								break;
							case "Размерность":
								album.isUseSize = true;
								buffItem.Size = Convert.ToString((worKsheeT.Cells[line, collCount] as Microsoft.Office.Interop.Excel.Range).Value);
								break;
							case "Кол-во":
								buffItem.Count = Convert.ToDouble((worKsheeT.Cells[line, collCount] as Microsoft.Office.Interop.Excel.Range).Value);
								break;
							case "Цена": 
								buffItem.Cost = Convert.ToDouble((worKsheeT.Cells[line, collCount] as Microsoft.Office.Interop.Excel.Range).Value);

								if (buffItem != null && buffItem.Count != 0 && buffItem.Cost != 0)
									album.Peoples[album.Peoples.Count - 1].Product.Add(buffItem);
								if (buffItem != null && buffItem.Count != 0 && buffItem.Cost == 0)
									throw new BadFormatException("Поле Цена незаполнено");
								buffItem = new SendProduct();
								break;
							default:
								throw new UnknownColumnException();
						}
						collCount++;

					}
					line++;
				}

				for (int i = 0; i < album.Peoples.Count; i++) {
					if (album.Peoples[i].Product.Count == 0) {
						album.Peoples.RemoveAt(i);
						i = i - 1;
					}

				}

			}
			catch (Exception ex) {
				throw ex;
			}
			finally {
				IntPtr handle = (IntPtr)excel.Hwnd;
				worKbooK.Close();
				excel.Quit();
				MyFunctions.TerminateProcess(handle, 1);
			}

			return album;
		}


		public static void deleteComments(int selectedIndex, VkNet.Utils.VkCollection<PhotoAlbum> AllAlbums, VkApi api, DateTime TimeFrom2, DateTime TimeTo2) {
			ulong offset = 0, count = 500;
            List<VkNet.Model.Attachments.Photo> AllPhoto = new List<VkNet.Model.Attachments.Photo>();
            try {
                while (true) {
                    var photos = api.Photo.Get(new PhotoGetParams { OwnerId = AllAlbums[selectedIndex].OwnerId, Extended = true, Count = count, Offset = offset, AlbumId = VkNet.Enums.SafetyEnums.PhotoAlbumType.Id(AllAlbums[selectedIndex].Id) });
                    for (int i = 0; i < photos.Count; i++) {
                        AllPhoto.Add(photos[i]);
                    }
                    if (AllPhoto.Count == AllAlbums[selectedIndex].Size)
                        break;
                    else if (AllPhoto.Count > AllAlbums[selectedIndex].Size) {
                        throw new OutOfBorderPhotoException();
                    }
                    else {
                        offset = offset + count;
                    }
                }
            }
            catch (Exception ex) {
                throw new ErrorFindPhotoException(AllAlbums[selectedIndex].Title, "ErrorFindPhotoException. " + ex.ToString());
            }

			for (int i = 0; i < AllPhoto.Count; i++) {
				if (AllPhoto[i].Comments.Count != 0) {
					count = 100;
					offset = 0;
					int delete_comment = 0;
					bool stop = true;
					try {
						while (true) {
							var Comments = api.Photo.GetComments(new PhotoGetCommentsParams { Sort = VkNet.Enums.SafetyEnums.CommentsSort.Desc, Count = count, Offset = offset, OwnerId = AllAlbums[selectedIndex].OwnerId, PhotoId = (ulong)AllPhoto[i].Id });
							for (int j = 0; j < Comments.Count; j++) {
								if (Comments[j].Date >= TimeFrom2 && Comments[j].Date <= TimeTo2) {
									api.Photo.DeleteComment(commentId: Convert.ToUInt64(Comments[j].Id), ownerId: -167157478);
									delete_comment++;
								}
								else {
									if (Comments[j].Date < TimeFrom2) {
										stop = false;
										break;
									}
								}
							}

							if (!stop || Comments.Count == 0)
								break;
							else {
								offset = offset + count;
							}
						}
					}
					catch (Exception ex){
						throw new ErrorFindCommentsException(i, AllAlbums[selectedIndex].Title, "ErrorFindCommentsException. " + ex.ToString());
					}
				}
			}
		}
	}
}
