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

namespace SoapOpera {
	/// <summary>
	/// Логика взаимодействия для MainWindow.xaml
	/// </summary>
	public partial class MainWindow : System.Windows.Window {
		VkApi api;
		const int groupId = 167157478;
		bool LoadCommplite = false;
		double AllCountPhoto = 0, CurrentCountPhoto = 0;
		VkNet.Utils.VkCollection<PhotoAlbum> AllAlbums;
		List<List<String>> AllNamePhoto;
		List<List<People>> peoples;
		DateTime TimeFrom, TimeTo;
		String ULogin = "79028669711";
		String ALogin = "79819776259";
		String UserPassword = "Folkner2017";
		String AntonPassword = "noskov-anton";

		private delegate void UpdateProgressBarDelegate(System.Windows.DependencyProperty dp, Object value);
		UpdateProgressBarDelegate updatePbDelegate;
		public MainWindow() {
			api = new VkApi();
			InitializeComponent();
			
			updatePbDelegate = new UpdateProgressBarDelegate(pbProgress.SetValue);

			TimeTo = DateTime.Now;
			tbDataTo.Text = TimeTo.Day + "." + TimeTo.Month + "." + TimeTo.Year;
			TimeTo2 = DateTime.Now;
			tbDataTo2.Text = TimeTo2.Day + "." + TimeTo2.Month + "." + TimeTo2.Year;
			try {
				LoadOptions();
				UserPassword = pbPassword.Password;
			}
			catch {
				LoadCommplite = true;
				lErrors.Content = "Не удалось загрузить параметры. Будут использованы параметры по умолчанию.";
			}

			try {
				LoadMessage();
			}
			catch {
				lErrors2.Document.Blocks.Clear();
				lErrors2.Document.Blocks.Add(new Paragraph(new Run("Не удалось загрузить сохранённые сообщения. Возможно файл с сообщиниями не существует.")));
				//lErrors2.Content = "Не удалось загрузить сохранённые сообщения. Возможно файл с сообщиниями не существует.";
			}

			try {
				api.Authorize(new ApiAuthParams {
					ApplicationId = 5539924,
					Login = ULogin,
					Password = UserPassword,
					Settings = Settings.All
				});
			
				uint countAlbums = (uint)api.Photo.GetAlbumsCount(groupId: groupId);
				AllAlbums = api.Photo.GetAlbums(new PhotoGetAlbumsParams { Count = countAlbums, OwnerId = -groupId }, skipAuthorization: true);
				for (int i = 0; i < countAlbums; i++) {
					cbAlbum.Items.Add("Альбом № " + (i + 1) + " имя: " + AllAlbums[i].Title + ". Количество фотографий: " + AllAlbums[i].Size);
					cbAlbum2.Items.Add("Альбом № " + (i + 1) + " имя: " + AllAlbums[i].Title + ". Количество фотографий: " + AllAlbums[i].Size);
				}
				cbAlbum.Items.Add("Все альбомы");
				cbAlbum.SelectedIndex = cbAlbum.Items.Count - 1;
				cbAlbum2.SelectedIndex = 0;
			}
			catch(Exception) {
				lErrors.Content = "Не удалось загрузить альбомы. Возможно нет подключения к интернету.";
			}
		}

		private void bPathFolder_Click(object sender, RoutedEventArgs e) {
			SaveFileDialog SaveFile = new SaveFileDialog();
			SaveFile.Filter = "xls files (*.xls)|*.xls";
			if (SaveFile.ShowDialog() != System.Windows.Forms.DialogResult.OK) {
				tbPahtFolder.Text = "";
			}
			else {
				tbPahtFolder.Text = SaveFile.FileName;
			}

		}

		private void bCreatFile_Click(object sender, RoutedEventArgs e) {
			if (ValidForm() == false) {
				return;
			}

			AllCountPhoto = CurrentCountPhoto = 0;
			List<string> NamesAlbums = new List<string>();
			peoples = new List<List<People>>();
			AllNamePhoto = new List<List<string>>();
			pbProgress.Value = 0;
			Dispatcher.Invoke(updatePbDelegate, System.Windows.Threading.DispatcherPriority.Background, new object[] { System.Windows.Controls.ProgressBar.ValueProperty, (double)0 });
			try {
				if (cbAlbum.SelectedIndex == cbAlbum.Items.Count - 1) {
					for (int i = 0; i < AllAlbums.Count; i++)
						AllCountPhoto = AllCountPhoto + (int)AllAlbums[i].Size;
					for (int i = 0; i < AllAlbums.Count; i++) {
						if(AllAlbums[i].Size != 0) {
							NamesAlbums.Add(AllAlbums[i].Title);
							peoples.Add(GetPeopleFromAlbum(i, api));
						}
					}
				}
				else {
					AllCountPhoto = (long)AllAlbums[cbAlbum.SelectedIndex].Size;
					NamesAlbums.Add(AllAlbums[cbAlbum.SelectedIndex].Title);
					peoples.Add(GetPeopleFromAlbum(cbAlbum.SelectedIndex, api));
				}
			}
			catch (OutOfBorderPhotoException) {
				lErrors.Content = "Выход за границу массива фотографий. Критическая ошибка.";
				return;
			}
			catch (OutOfBorderCommentsException) {
				lErrors.Content = "Выход за границу массива комментариев. Критическая ошибка.";
				return;
			}
			catch (ErrorFindPhotoException exc) {
				lErrors.Content = "Не удалось получить все фотографии из альбома " + exc.TitleAlbum + ". Возможно отсутствует подключение к интернету.";
				return;
			}
			catch (ErrorFindCommentsException exc) {
				lErrors.Content = "Не удалось получить коментарии к фотографии № " + exc.NumberPhoto + " из альбома " + exc.TitleAlbum + ". Возможно отсутствует подключение к интернету.";
				return;
			}
			catch (ErrorFindNameUserException exp) {
				lErrors.Content = "Не удалось получить имена пользователей в альбоме " + exp.TitleAlbum + " . Возможно отсутствует подключение к интернету.";
				return;
			}
			catch (Exception) {
				lErrors.Content = "Критическая ошибка.";
				return;
			}

			try {
				if (cbUseColColum.IsChecked.Value) {
					MyFunctions.WriteExcelFileWithCollum(peoples, NamesAlbums, cbExeptEmptyPhoto.IsChecked.Value, AllNamePhoto, cbFormat.SelectedIndex, tbPahtFolder.Text);
				}
				else {
					MyFunctions.WriteExcelFile(peoples, NamesAlbums, cbExeptEmptyPhoto.IsChecked.Value, AllNamePhoto, cbFormat.SelectedIndex, tbPahtFolder.Text);
				}
			}
			catch (NotCommentException){
				lErrors.Content = "Файл не был создан, так как не нашлось комментариев для записи.";
				return;
			}
			catch (Exception) {
				lErrors.Content = "Критическая ошибка при сохранении файла.";
				return;
			}

			lErrors.Content = "Файл сформирован.";
		}

		List<People> GetPeopleFromAlbum(int SelectedIndex, VkApi api) {
			ulong offset = 0, count = 500;
			List<VkNet.Model.Attachments.Photo> AllPhoto = new List<VkNet.Model.Attachments.Photo>();
			List<List<VkNet.Model.Comment>> AllComents = new List<List<VkNet.Model.Comment>>();
			AllNamePhoto.Add(new List<string>());


			try {
				while (true) {
					var photos = api.Photo.Get(new PhotoGetParams { OwnerId = AllAlbums[SelectedIndex].OwnerId, Extended = true, Count = count, Offset = offset, AlbumId = VkNet.Enums.SafetyEnums.PhotoAlbumType.Id(AllAlbums[SelectedIndex].Id) });
					for (int i = 0; i < photos.Count; i++) {
						AllNamePhoto[AllNamePhoto.Count - 1].Add(photos[i].Text.Split('.')[0].Split('\n')[0]);
						AllPhoto.Add(photos[i]);
					}
					if (AllPhoto.Count == AllAlbums[SelectedIndex].Size)
						break;
					else if (AllPhoto.Count > AllAlbums[SelectedIndex].Size) {
						throw new OutOfBorderPhotoException();
					}
					else {
						offset = offset + count;
					}
				}
			}
			catch {
				throw new ErrorFindPhotoException(AllAlbums[SelectedIndex].Title);
			}

			for (int i = 0; i < AllPhoto.Count; i++) {
				CurrentCountPhoto = CurrentCountPhoto + 1;
				Dispatcher.Invoke(updatePbDelegate, System.Windows.Threading.DispatcherPriority.Background, new object[] { System.Windows.Controls.ProgressBar.ValueProperty, (double)CurrentCountPhoto * (double)100 / (double)AllCountPhoto });
				AllComents.Add(new List<VkNet.Model.Comment>());
				if (AllPhoto[i].Comments.Count != 0) {
					count = 100;
					offset = 0;
					bool stop = true;
					try {
						while (true) {
							var Comments = api.Photo.GetComments(new PhotoGetCommentsParams { Sort = VkNet.Enums.SafetyEnums.CommentsSort.Desc, Count = count, Offset = offset, OwnerId = AllAlbums[SelectedIndex].OwnerId, PhotoId = (ulong)AllPhoto[i].Id });
							for (int j = 0; j < Comments.Count; j++) {
								if (Comments[j].Date >= TimeFrom && Comments[j].Date <= TimeTo) {
									AllComents[i].Add(Comments[j]);
								}
								else {
									if (Comments[j].Date < TimeFrom) {
										stop = false;
										break;
									}
								}
							}
							if (AllComents[i].Count == AllPhoto[i].Comments.Count || !stop || Comments.Count == 0)
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
						throw new ErrorFindCommentsException(i, AllAlbums[SelectedIndex].Title);
					}
				}
			}
			List<People> peoples = new List<People>();
			offset = 0;
			count = 500;
			List<List<long>> Users = new List<List<long>>();
			Users.Add(new List<long>());
			for (int i = 0; i < AllComents.Count; i++) {
				for (int j = 0; j < AllComents[i].Count; j++) {
					People buff;
					if (cbExeptGroupId.IsChecked == false || AllComents[i][j].FromId != -groupId) {
						if ((buff = peoples.Find(x => x.Id == AllComents[i][j].FromId)) != null) {
							buff.CommentsText[i].Add(AllComents[i][j].Text);
							buff.CommentsData[i].Add(AllComents[i][j].Date.Value);
						}
						else {
							peoples.Add(new People(AllComents.Count, "", AllComents[i][j].FromId));
							if (Users[(int)offset].Count == (int)count) {
								Users.Add(new List<long>());
								offset++;
							}
							Users[(int)offset].Add(AllComents[i][j].FromId);

							peoples[peoples.Count - 1].CommentsText[i].Add(AllComents[i][j].Text);
							peoples[peoples.Count - 1].CommentsData[i].Add(AllComents[i][j].Date.Value);
						}
					}
				}
			}
			try {
				int k = 0;
				for (int i = 0; i < Users.Count; i++) {
					if (Users[i].Count == 0)
						continue;
					var UsersVk = api.Users.Get(userIds: Users[i]);
					for (int j = 0; j < UsersVk.Count; j++) {
						if (peoples[k].Id == -groupId) {
							var GroupVk = api.Groups.GetById(new string[] { (groupId).ToString() }, (groupId).ToString(), GroupsFields.CityId);
							peoples[k].Name = GroupVk[0].Name;
							j--;
						}
						else
							peoples[k].Name = UsersVk[j].FirstName + " " + UsersVk[j].LastName;
						k++;
					}
				}
			}
			catch (Exception) {
				throw new ErrorFindNameUserException(AllAlbums[SelectedIndex].Title);
			}
			return peoples;
		}

		private void cbExeptGroupId_Checked(object sender, RoutedEventArgs e) {
			try {
				if (LoadCommplite)
					SaveOptions();
			}
			catch {
				lErrors.Content = "Не удалось сохранить параметры.";
			}
		}

		private void cbExeptEmptyPhoto_Checked(object sender, RoutedEventArgs e) {
			try {
				if (LoadCommplite)
					SaveOptions();
			}
			catch {
				lErrors.Content = "Не удалось сохранить параметры.";
			}
		}

		void SaveOptions() {
			XmlDocument doc = new XmlDocument();

			doc.AppendChild(doc.CreateElement("Options"));

			XmlElement CounterXml = doc.CreateElement("CreateDocument");
			CounterXml.SetAttribute("Exept_Group_Id", cbExeptGroupId.IsChecked.ToString());
			CounterXml.SetAttribute("Exept_Empty_Photo", cbExeptEmptyPhoto.IsChecked.ToString());
			doc.ChildNodes[0].AppendChild(CounterXml);

			CounterXml = doc.CreateElement("UserPassword");
			CounterXml.SetAttribute("Password", pbPassword.Password);
			doc.ChildNodes[0].AppendChild(CounterXml);

			doc.Save(Environment.CurrentDirectory + "\\Options.xml");
		}

		private void pbPassword_PasswordChanged(object sender, RoutedEventArgs e) {
			try {
				if (LoadCommplite)
					SaveOptions();
			}
			catch {
				lErrors.Content = "Не удалось сохранить параметры.";
			}
		}

		void LoadOptions() {
			XmlDocument doc = new XmlDocument();

			doc.Load(Environment.CurrentDirectory + "\\Options.xml");
			XmlNode CounterXml = doc.ChildNodes[0].ChildNodes[0];

			cbExeptGroupId.IsChecked = (CounterXml.Attributes[0].Value == "True");
			cbExeptEmptyPhoto.IsChecked = (CounterXml.Attributes[1].Value == "True");

			CounterXml = doc.ChildNodes[0].ChildNodes[1];

			pbPassword.Password = CounterXml.Attributes[0].Value;

			LoadCommplite = true;
		}

		bool ValidForm() {
			lErrors.Content = "";
			if (tbPahtFolder.Text.Length == 0) {
				lErrors.Content = "Не выбран файл для сохранения.";
				return false;
			}

			if (cbAlbum.SelectedIndex == -1) {
				lErrors.Content = "Не выбран альбом.";
				return false;
			}

			try {
				string[] splitString = tbDataFrom.Text.Split(new Char[] { '.' }, StringSplitOptions.RemoveEmptyEntries);
				if (splitString.Length == 3)
					TimeFrom = new DateTime(Convert.ToInt32(splitString[2]), Convert.ToInt32(splitString[1]), Convert.ToInt32(splitString[0]));
				else
					throw new Exception();
			}
			catch {
				lErrors.Content = "Не правильно введена дата, с которой начинается поиск. Формат dd.mm.yyyy";
				return false;
			}

			try {
				string[] splitString = tbDataTo.Text.Split(new Char[] { '.' }, StringSplitOptions.RemoveEmptyEntries);
				if (splitString.Length == 3)
					TimeTo = new DateTime(Convert.ToInt32(splitString[2]), Convert.ToInt32(splitString[1]), Convert.ToInt32(splitString[0]),23,59,59);
				else
					throw new Exception();
			}
			catch {
				lErrors.Content = "Не правильно введена дата, до которой идёт поиск. Формат dd.mm.yyyy";
				return false;
			}

			if(TimeFrom > TimeTo) {
				lErrors.Content = "Дата, до которой идёт поиск, должна быть больше даты начала поиска.";
				return false;
			}

			if (cbAlbum.SelectedIndex  != cbAlbum.Items.Count - 1 && AllAlbums[cbAlbum.SelectedIndex].Size == 0) {
				lErrors.Content = "Не возможно сформировать файл. В альбоме нет фотографий.";
				return false;
			}
			Microsoft.Office.Interop.Excel.Application excel = null;
			try {
				excel = new Microsoft.Office.Interop.Excel.Application();
				excel.Visible = false;
				excel.DisplayAlerts = false;
				Workbook worKbooK = excel.Workbooks.Add(Type.Missing);

				if (cbFormat.SelectedIndex == 0 || cbFormat.SelectedIndex == -1)
					worKbooK.SaveAs(tbPahtFolder.Text, XlFileFormat.xlExcel3);
				else
					worKbooK.SaveAs(tbPahtFolder.Text);

				worKbooK.Close();
			}
			catch (Exception exc) {
				lErrors.Content = "Не удаётся создать или пересохранить файл. Выберите другое имя.";
				return false;
			}
			finally {
				IntPtr handle = (IntPtr)excel.Hwnd;
				excel.Quit();
				MyFunctions.TerminateProcess(handle, 1);
			}

			return true;
		}
	}
}
