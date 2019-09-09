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

		private void bResiveMessage_Click(object sender, RoutedEventArgs e) {
			if (!ValidForm2())
				return;

			SendAlbum album;
			//Формирование списка людей
			try {
				lErrors2.Document.Blocks.Clear();
				album = MyFunctions.GetSendPeoples(tbPathFolder2.Text);
			}
			catch (BadFormatException ex) {
				lErrors2.Document.Blocks.Add(new Paragraph(new Run(ex.Error)));
				return;
			}
			catch (UnknownColumnException ex) {
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Найден неизвестный столбец")));
				return;
			}
			catch (Exception ex) {
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось прочесть файл. Возможно не хватает некоторых полей.")));
				return;
			}

			//Авторатизация в vk
			var api2 = new VkApi();

			if (!api2.IsAuthorized) {
				try {
					api2.Authorize(new ApiAuthParams {
						ApplicationId = 5539924,
                        //Login = ULogin,
                        //Login = ALogin,

                        //Password = UserPassword,
                        //Password = AntonPassword,

                        AccessToken = "c2d23fe54097f46ea81f8103efd840bc5ca545408c14ee8873c9ee83711d123ca2d5a049f57b774028dfd",
                        Settings = Settings.All
					});
				}
				catch (Exception ex) {
					lErrors2.Document.Blocks.Clear();
					lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось авторитизироваться в VK. Возможно нет подключения к интеренту. " + ex.ToString())));
					//lErrors2.Content = "Не удалось авторитизироваться в VK. Возможно нет подключения к интеренту.";
					return;
				}
			}
			//lErrors2.Content = "";
			//Формирование и отправка сообщений
			for (int i = 0; i < album.Peoples.Count; i++) {
				string Message = new TextRange(rtbTextArea.Document.ContentStart, rtbTextArea.Document.ContentEnd).Text;
				Message = Message.Replace("{имя}", album.Peoples[i].Name);

				if (rbChouse1.IsChecked == true) {
					Message = Message.Replace("{стоимость}", album.Peoples[i].AllCost.ToString());
					Message = Message.Replace("{альбом}", album.Name);
					string listProduct = "";
					for (int j = 0; j < album.Peoples[i].Product.Count; j++) {
						listProduct = listProduct + (j + 1) + ". " + album.Peoples[i].Product[j].Name + ((album.Peoples[i].Product[j].Size == "") ? "" : " " + album.Peoples[i].Product[j].Size + " ") + ": " + album.Peoples[i].Product[j].Count + " шт * " + album.Peoples[i].Product[j].Cost + " руб = " + album.Peoples[i].Product[j].Count * album.Peoples[i].Product[j].Cost + " руб.\n";
					}
					Message = Message.Replace("{список}", listProduct);
				}


				long? captcha_sid = null;
				string captcha_key = null;
                Random prob1 = new Random();

                for (int j = 0; j < 10; j++) {
					try {
						api2.Messages.Send(new MessagesSendParams { UserId = album.Peoples[i].Id, Message = Message, CaptchaSid = captcha_sid, CaptchaKey = captcha_key, RandomId = prob1.Next() });
						break;
					}
					catch (CaptchaNeededException ex) {
						CaptchaForm prob = new CaptchaForm(ex.Img, "Сообщение № " + i);
						prob.ShowDialog();
						captcha_sid = ex.Sid;
						captcha_key = prob.CaptchaKey;
						if (j == 9)
							throw new Exception();

					}
					catch (Exception exc) {
						lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось отправить сообщение покупателю " + album.Peoples[i].Name + ".")));
						if(exc.Message == "Can't send messages for users without permission") {
                            lErrors2.Document.Blocks.Add(new Paragraph(new Run("У пользователя запрещены сообщения от группы.")));
                        }
                        lErrors2.Document.Blocks.Add(new Paragraph(new Run(exc.ToString())));
                        //lErrors2.Content = lErrors2.Content + "Не удалось отправить сообщение покупателю " + peoples[i].Name + ". " + exc.ToString();
                        break;
					}
				}
			}


			if (new TextRange(lErrors2.Document.ContentStart, lErrors2.Document.ContentEnd).Text == "") {
				lErrors2.Document.Blocks.Clear();
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Сообщения отправлены.")));
				//lErrors2.Content = "Сообщения отправлены.";
			}

		}

		private void bPathFolder2_Click(object sender, RoutedEventArgs e) {
			SaveFileDialog SaveFile = new SaveFileDialog();
			SaveFile.Filter = "xls files (*.xls)|*.xls";
			tbPathFolder2.Text = (SaveFile.ShowDialog() == System.Windows.Forms.DialogResult.OK) ? SaveFile.FileName : "";
		}

		private void bSaveMassage_Click(object sender, RoutedEventArgs e) {
			lErrors2.Document.Blocks.Clear();
			//lErrors2.Content = "";
			try {
				SaveMessage();
			}
			catch {
				lErrors2.Document.Blocks.Clear();
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось сохранить сообщение.")));
				//lErrors2.Content = "Не удалось сохранить сообщение.";
			}
		}

		private void rbChouse1_Checked(object sender, RoutedEventArgs e) {
			if (rtbTextArea != null) {
				try {
					LoadMessage();
				}
				catch {
					lErrors2.Document.Blocks.Clear();
					lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось загрузить сохранённые сообщения. Воз можно файл с сообщиниями отсуцтвует.")));
					//lErrors2.Content = "Не удалось загрузить сохранённые сообщения. Воз можно файл с сообщиниями отсуцтвует.";
				}
			}
		}

		private void rbChouse2_Checked(object sender, RoutedEventArgs e) {
			if (rtbTextArea != null) {
				try {
					LoadMessage();
				}
				catch {
					lErrors2.Document.Blocks.Clear();
					lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось загрузить сохранённые сообщения. Возможно файл с сообщиниями отсуцтвует.")));
					//lErrors2.Content = "Не удалось загрузить сохранённые сообщения. Возможно файл с сообщиниями отсуцтвует.";
				}
			}
		}

		void SaveMessage() {
			XmlDocument doc = new XmlDocument();
			string CountMessage = "", OtherMessage = "";

			if (rbChouse1.IsChecked == true)
				CountMessage = new TextRange(rtbTextArea.Document.ContentStart, rtbTextArea.Document.ContentEnd).Text;
			else
				OtherMessage = new TextRange(rtbTextArea.Document.ContentStart, rtbTextArea.Document.ContentEnd).Text;

			try {
				doc.Load(Environment.CurrentDirectory + "\\Messages.xml");
				if(rbChouse1.IsChecked == true)
					OtherMessage = doc.ChildNodes[0].ChildNodes[1].Attributes[0].Value;
				else
					CountMessage = doc.ChildNodes[0].ChildNodes[0].Attributes[0].Value;

				doc.ChildNodes[0].ChildNodes[0].Attributes[0].Value = CountMessage;
				doc.ChildNodes[0].ChildNodes[1].Attributes[0].Value = OtherMessage;
				doc.Save(Environment.CurrentDirectory + "\\Messages.xml");
				return;
			}
			catch {
				doc.AppendChild(doc.CreateElement("Messages"));

				XmlElement CounterXml = doc.CreateElement("Message");
				CounterXml.SetAttribute("CountMessage", CountMessage);
				doc.ChildNodes[0].AppendChild(CounterXml);

				CounterXml = doc.CreateElement("Message");
				CounterXml.SetAttribute("OtherMessage", OtherMessage);
				doc.ChildNodes[0].AppendChild(CounterXml);

				doc.Save(Environment.CurrentDirectory + "\\Messages.xml");
			}
		}



		void LoadMessage() {
			XmlDocument doc = new XmlDocument();
			doc.Load(Environment.CurrentDirectory + "\\Messages.xml");
			string buff; 
			if (rbChouse1.IsChecked == true) {
				buff = doc.ChildNodes[0].ChildNodes[0].Attributes[0].Value;

			}
			else {
				buff = doc.ChildNodes[0].ChildNodes[1].Attributes[0].Value;
			}
			rtbTextArea.Document.Blocks.Clear();
			rtbTextArea.Document.Blocks.Add(new Paragraph(new Run(buff)));
		}

		bool ValidForm2() {
			lErrors2.Document.Blocks.Clear();
			//lErrors2.Content = "";
			if (tbPathFolder2.Text.Length == 0) {
				lErrors2.Document.Blocks.Clear();
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не выбран файл для отправки.")));
				//lErrors2.Content = "Не выбран файл для сохранения.";
				return false;
			}

			if(!File.Exists(tbPathFolder2.Text)) {
				lErrors2.Document.Blocks.Clear();
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Файл не существует.Выберите другой файл.")));
				//lErrors2.Content = "Файл не существует.Выберите другой файл.";
				return false;
			}
			string prob = new TextRange(rtbTextArea.Document.ContentStart, rtbTextArea.Document.ContentEnd).Text.Replace("\n","").Replace("\r","").Replace(" ", "");
			if (prob.Length == 0) {
				lErrors2.Document.Blocks.Clear();
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Нет текста сообщения.")));
				//lErrors2.Content = "Нет текста сообщения.";
				return false;
			}

			return true;
		}

	}
}


/*
 * private void bResiveMessage_Click(object sender, RoutedEventArgs e) {
			if (!ValidForm2())
				return;

			SendAlbum p = GetSendPeoples(tbPathFolder2.Text);
			List<People> peoples = new List<People>();
			List<List<PhotoInfo>> photoInfo;
			String AlbumName = "";
			bool provColl = false;
			Microsoft.Office.Interop.Excel.Application excel = new Microsoft.Office.Interop.Excel.Application();
			Workbook worKbooK = null;
			excel.Visible = false;
			excel.DisplayAlerts = false;
			try {		
				worKbooK = excel.Workbooks.Open(tbPathFolder2.Text);
				Worksheet worKsheeT = (Worksheet)worKbooK.ActiveSheet;
				int countPhoto = 4;
				int line = 1;

				AlbumName = (string)(worKsheeT.Cells[line, 1] as Microsoft.Office.Interop.Excel.Range).Value;


				line = 3;
				if ((string)(worKsheeT.Cells[line, 5] as Microsoft.Office.Interop.Excel.Range).Value == "Размерность") {
					provColl = true;
				}

				line = 2;
				while (true) {
					string buff = (string)(worKsheeT.Cells[line, countPhoto] as Microsoft.Office.Interop.Excel.Range).Value;
					if (buff == "Сумма к оплате")
						break;
					else
						if (buff == null && (string)(worKsheeT.Cells[line, countPhoto + 1] as Microsoft.Office.Interop.Excel.Range).Value == null && (string)(worKsheeT.Cells[line, countPhoto + 2] as Microsoft.Office.Interop.Excel.Range).Value == null)
							throw new BadFormatException("Не найден заголовок \"Сумма к оплате\".");
						countPhoto++;
				}


				line = 4;
				photoInfo = new List<List<PhotoInfo>>();
				while (true) {
					string name = null;
					long id = 0;
					double Buffcount = 0; 

					try {
						name = (string)(worKsheeT.Cells[line, 2] as Microsoft.Office.Interop.Excel.Range).Value;
					}
					catch (Exception ex) {
						throw new BadFormatException("Неправельное имя пользователя. Строка: " + line + ".");
					}

					try {
						id = Convert.ToInt64((worKsheeT.Cells[line, 3] as Microsoft.Office.Interop.Excel.Range).Value);
						if (id == 0 && name != null)
							throw new Exception();
					}
					catch {
						throw new BadFormatException("Неправельное id пользователя. Строка: " + line + ".");
					}

					if (rbChouse1.IsChecked == true) {

						photoInfo.Add(new List<PhotoInfo>());
						for (int i = 5;i < countPhoto;i = i + 3) {
							try {
								double count = Convert.ToDouble((worKsheeT.Cells[line, i] as Microsoft.Office.Interop.Excel.Range).Value);
								double price = Convert.ToDouble((worKsheeT.Cells[line, i + 1] as Microsoft.Office.Interop.Excel.Range).Value);
								if (count != 0 && price != 0)  {
									photoInfo[photoInfo.Count - 1].Add(new PhotoInfo(count,price, (string)(worKsheeT.Cells[2, i - 1] as Microsoft.Office.Interop.Excel.Range).Value));
								}
							}
							catch {}
						}


						try {
							Buffcount = (worKsheeT.Cells[line, countPhoto] as Microsoft.Office.Interop.Excel.Range).Value;
							if (Buffcount == 0)
								throw new Exception();
						}
						catch (Exception ex) {
							throw new BadFormatException("Неправельное поле \"Сумма к оплате\". Строка: " + line + ". Поле \"Сумма к оплате\" ровняется нулю.");
						}
					}
					else {
						Buffcount = 0;
					}

					if (rbChouse2.IsChecked == true || photoInfo[line-4].Count != 0) {
						peoples.Add(new People(0, name, id));
						peoples[peoples.Count - 1].Count = Buffcount;
						line++;
					}
					else {
						throw new BadFormatException("Не одного товара у покупателя. Строка " + line + ".");
					}

					if ((string)(worKsheeT.Cells[line, 2] as Microsoft.Office.Interop.Excel.Range).Value == null) {
						break;
					}
				}
			}
			catch(BadFormatException E) {
				lErrors2.Document.Blocks.Clear();
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось прочесть файл. " + E.Error)));
				//lErrors2.Content = "Не удалось прочесть файл. " + E.Error;
				return;
			}
			catch (Exception) {
				lErrors2.Document.Blocks.Clear();
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось прочесть файл. Возможно не хватает некоторых полей.")));
				//lErrors2.Content = "Не удалось прочесть файл. Возможно не хватает некоторых полей.";
				return;
			}
			finally {
				IntPtr handle = (IntPtr)excel.Hwnd;
				worKbooK.Close();
				excel.Quit();
				MyFunctions.TerminateProcess(handle,1);
			}

			var api2 = new VkApi();

			if (!api2.IsAuthorized) {
				try {
					api2.Authorize(new ApiAuthParams {
						ApplicationId = 5539924,
						Login = ALogin,
						Password = AntonPassword,
						Settings = Settings.All
					});
				}
				catch {
					lErrors2.Document.Blocks.Clear();
					lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось авторитизироваться в VK. Возможно нет подключения к интеренту.")));
					//lErrors2.Content = "Не удалось авторитизироваться в VK. Возможно нет подключения к интеренту.";
					return;
				}
			}
			lErrors2.Document.Blocks.Clear();
			//lErrors2.Content = "";
			for (int i = 0; i < peoples.Count;i++) {
				string Message = new TextRange(rtbTextArea.Document.ContentStart, rtbTextArea.Document.ContentEnd).Text;
				Message = Message.Replace("{имя}", peoples[i].Name);

				if (rbChouse1.IsChecked == true) {
					Message = Message.Replace("{стоимость}", peoples[i].Count.ToString());
					Message = Message.Replace("{альбом}", AlbumName);
					string listProduct = "";
					for (int j = 0; j < photoInfo[i].Count;j++) {
						listProduct = listProduct + (j + 1) + ". " + photoInfo[i][j].PhotoName + ": " + photoInfo[i][j].Count + " шт * " + photoInfo[i][j].Price + " руб = " +photoInfo[i][j].Count * photoInfo[i][j].Price + " руб.\n";
					}
					Message = Message.Replace("{список}",listProduct);
				}
				long? captcha_sid = null;
				string captcha_key = null;
				try {
					for (int j = 0; j < 10; j++) {
						try {
							api2.Messages.Send(new MessagesSendParams { UserId = peoples[i].Id, Message = Message, CaptchaSid = captcha_sid, CaptchaKey = captcha_key });
							break;
						}
						catch (CaptchaNeededException ex) {
							CaptchaForm prob = new CaptchaForm(ex.Img,"Сообщение № " + i);
							prob.ShowDialog();
							captcha_sid = ex.Sid;
							captcha_key = prob.CaptchaKey;
							if (j == 9)
								throw new Exception();

						}
						catch (Exception exc) {
							lErrors2.Document.Blocks.Add(new Paragraph(new Run(new TextRange(lErrors2.Document.ContentStart, lErrors2.Document.ContentEnd).Text + "Не удалось отправить сообщение покупателю " + peoples[i].Name + ". \n" + exc.ToString())));
							//lErrors2.Content = lErrors2.Content + "Не удалось отправить сообщение покупателю " + peoples[i].Name + ". " + exc.ToString();
							break;
						}
					}
				}
				catch (Exception E) {
					lErrors2.Document.Blocks.Add(new Paragraph(new Run(new TextRange(lErrors2.Document.ContentStart, lErrors2.Document.ContentEnd).Text + "Не удалось отправить сообщение покупателю " + peoples[i].Name + ". \n")));
				}
			}
			if (new TextRange(lErrors2.Document.ContentStart, lErrors2.Document.ContentEnd).Text == "") {
				lErrors2.Document.Blocks.Clear();
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Сообщения отправлены.")));
				//lErrors2.Content = "Сообщения отправлены.";
			}
		}

	*/