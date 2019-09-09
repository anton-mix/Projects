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
using System.Runtime.InteropServices;

using Microsoft.Office.Interop.Excel;
using System.Xml;

using VkNet;
using VkNet.Exception;
using VkNet.Enums.Filters;
using VkNet.Infrastructure;
using VkNet.Model;
using VkNet.Model.RequestParams;

namespace SoapOpera {
	public partial class MainWindow : System.Windows.Window {


		private void CreatePrintFile_Click(object sender, RoutedEventArgs e) {
			if (!ValidPrintForm())
				return;
			SendAlbum album;
			try {
				lErrors3.Document.Blocks.Clear();
				album = MyFunctions.GetSendPeoples(tbPathFolderPrint1.Text);
			}
			catch (BadFormatException ex) {
				lErrors3.Document.Blocks.Add(new Paragraph(new Run(ex.Error)));
				return;
			}
			catch (UnknownColumnException ex) {
				lErrors3.Document.Blocks.Add(new Paragraph(new Run("Найден неизвестный столбец")));
				return;
			}
			catch (Exception ex) {
				lErrors3.Document.Blocks.Add(new Paragraph(new Run("Не удалось прочесть файл. Возможно не хватает некоторых полей. " + ex.ToString())));
				return;
			}

			try {
				lErrors3.Document.Blocks.Clear();
				CreatePrintFiles(album, tbPathFolderPrint2.Text, cbFormat.SelectedIndex);
				lErrors3.Document.Blocks.Add(new Paragraph(new Run("Файл сформирован.")));
			}
			catch (BadFormatSafeException ex) {
				lErrors3.Document.Blocks.Add(new Paragraph(new Run("Невозможно сохранить. Возможно файл открыт или неправельно выбран формат Excel.")));
			}
			catch (Exception ex) {
				lErrors3.Document.Blocks.Add(new Paragraph(new Run("Произошла ошибка при создании файла. " + ex.ToString())));
			}

		}

		private void CreatePrintFiles(SendAlbum album, String NameFile, int Format) {
			int line = 1;

			Microsoft.Office.Interop.Excel.Application excel = new Microsoft.Office.Interop.Excel.Application();
			excel.Visible = false;
			excel.DisplayAlerts = false;
			Workbook worKbooK = excel.Workbooks.Add(Type.Missing);
			Worksheet worKsheeT = (Worksheet)worKbooK.ActiveSheet;
			worKsheeT.Name = "SoapOperaPrint";
			try {
				line = 1;
				worKsheeT.Cells[line, 1] = album.Name;

				line = 2;
				worKsheeT.Cells[line, 1] = "№";
				worKsheeT.Cells[line, 2] = "Участник";

				if (album.isUseSize)
					for (int i = 3; i < album.Photos.Count * 2 + 3; i = i + 2) {
						worKsheeT.Cells[line, i] = album.Photos[(i - 3) / 2];
						worKsheeT.Cells[line, i].Orientation = 90;
					}
				else
					for (int i = 3; i < album.Photos.Count + 3; i++) {
						worKsheeT.Cells[line, i] = album.Photos[i - 3];
						worKsheeT.Cells[line, i].Orientation = 90;
					}


				line = 3;
				for (int i = 0; i < album.Peoples.Count; i++) {
					worKsheeT.Cells[line, 1] = (i + 1);
					worKsheeT.Cells[line, 2] = album.Peoples[i].Name;
					for (int j = 0; j < album.Peoples[i].Product.Count; j++) {
						for (int k = 0; k < album.Photos.Count; k++) {
							if (album.Peoples[i].Product[j].Name == album.Photos[k]) {
								if (album.isUseSize) {
									worKsheeT.Cells[line, 3 + k * 2] = album.Peoples[i].Product[j].Size;
									worKsheeT.Cells[line, 3 + k * 2 + 1] = album.Peoples[i].Product[j].Count;
								}
								else {
									worKsheeT.Cells[line, 3 + k] = album.Peoples[i].Product[j].Count;
								}
							}
						}
					}
					line++;
				}

				int countCollum = 0;
				if(album.isUseSize) {
					countCollum = 2 + album.Photos.Count * 2;
				}
				else {
					countCollum = 2 + album.Photos.Count;
				}

				Microsoft.Office.Interop.Excel.Borders border = worKsheeT.Range[worKsheeT.Cells[2, 1], worKsheeT.Cells[line - 1, countCollum]].Borders;
				border.LineStyle = Microsoft.Office.Interop.Excel.XlLineStyle.xlContinuous;
				border.Weight = 2d;


				Microsoft.Office.Interop.Excel.Range usedrange = worKsheeT.Range[worKsheeT.Cells[1, 1], worKsheeT.Cells[line-1, countCollum]];
				usedrange.Rows.AutoFit();
				usedrange.Columns.AutoFit();

				try {
					if (Format <= 0)
						worKbooK.SaveAs(NameFile, XlFileFormat.xlExcel3);
					else
						worKbooK.SaveAs(NameFile);
				}
				catch (Exception ex2) {
					throw new BadFormatSafeException();
				}
				worKbooK.Close();
			}
			catch (Exception ex) {
				throw ex;
			}
			finally {
				IntPtr handle = (IntPtr)excel.Hwnd;
				excel.Quit();
				MyFunctions.TerminateProcess(handle, 1);
			}
		}

		private void bPathFolderPrint1_Click(object sender, RoutedEventArgs e) {
			SaveFileDialog SaveFile = new SaveFileDialog();
			SaveFile.Filter = "xls files (*.xls)|*.xls";
			tbPathFolderPrint1.Text = (SaveFile.ShowDialog() == System.Windows.Forms.DialogResult.OK) ? SaveFile.FileName : "";
		}

		private void bPathFolderPrint2_Click(object sender, RoutedEventArgs e) {
			SaveFileDialog SaveFile = new SaveFileDialog();
			SaveFile.Filter = "xls files (*.xls)|*.xls";
			tbPathFolderPrint2.Text = (SaveFile.ShowDialog() == System.Windows.Forms.DialogResult.OK) ? SaveFile.FileName : "";
		}

		private bool ValidPrintForm() {
			lErrors3.Document.Blocks.Clear();
			if (tbPathFolderPrint1.Text == "") {
				lErrors3.Document.Blocks.Add(new Paragraph(new Run("Не выбран файл с таблицей.")));
				return false;
			}
			if (tbPathFolderPrint2.Text == "") {
				lErrors3.Document.Blocks.Add(new Paragraph(new Run("Не выбран файл для печати.")));
				return false;
			}
			return true;
		}

	}
}
